/// myco_downloader
///
/// Downloads Mycobacterium gene sequences from NCBI using the Entrez
/// E-utilities API (esearch + efetch).
///
/// Supported targets: 16s, hsp65, rpob, erm41
///
/// Usage:
///   cargo run -- --gene 16s  --output myco_16s.fasta
///   cargo run -- --gene hsp65
///   cargo run -- --gene rpob --max 500 --email your@email.com
use anyhow::{Context, Result};
use clap::Parser;
use indicatif::{ProgressBar, ProgressStyle};
use std::fs::OpenOptions;
use std::io::Write;
use std::thread;
use std::time::Duration;

// ── CLI ──────────────────────────────────────────────────────────────────────

#[derive(Parser, Debug)]
#[command(
    name = "myco_downloader",
    about = "Download Mycobacterium gene sequences from NCBI Entrez"
)]
struct Args {
    /// Gene target: 16s | hsp65 | rpob | erm41
    #[arg(short, long, default_value = "16s")]
    gene: String,

    /// Output FASTA file path (defaults to myco_<gene>.fasta)
    #[arg(short, long)]
    output: Option<String>,

    /// Maximum sequences to download (0 = all found)
    #[arg(short, long, default_value_t = 0)]
    max: usize,

    /// Your email (required by NCBI for API usage)
    #[arg(short, long, default_value = "researcher@example.com")]
    email: String,

    /// NCBI API key (optional — increases rate limit from 3 to 10 req/s)
    #[arg(short, long)]
    api_key: Option<String>,

    /// Batch size per efetch call (max 500 recommended)
    #[arg(short, long, default_value_t = 200)]
    batch: usize,

    /// Print summary stats only, do not download sequences
    #[arg(long)]
    count_only: bool,
}

// ── Gene targets ─────────────────────────────────────────────────────────────

struct Target {
    label: &'static str,
    query: &'static str,
    default_output: &'static str,
}

fn resolve_target(gene: &str) -> Result<Target> {
    match gene.to_lowercase().as_str() {
        "rrs" => Ok(Target {
            label: "16S rRNA",
            query: "Mycobacteriaceae[Organism] AND (16S[Title] OR rrs[Gene Name]) \
                    AND 400:1800[SLEN] AND biomol_rrna[PROP]",
            default_output: "myco_rrs.fasta",
        }),
        "hsp65" => Ok(Target {
            label: "hsp65",
            query: "Mycobacteriaceae[Organism] AND (hsp65[Gene Name] OR groEL2[Gene Name]) \
                    AND 300:700[SLEN]",
            default_output: "myco_hsp65.fasta",
        }),
        "rpob" => Ok(Target {
            label: "rpoB",
            query: "Mycobacteriaceae[Organism] AND rpoB[Gene Name] AND 300:800[SLEN]",
            default_output: "myco_rpob.fasta",
        }),
        "erm41" => Ok(Target {
            label: "erm(41)",
            query: "Mycobacteriaceae[Organism] AND erm(41)[Gene Name] AND 300:1000[SLEN]",
            default_output: "myco_erm41.fasta",
        }),
        "rrl" => Ok(Target {
            label: "23S rRNA (rrl)",
            query: "Mycobacteriaceae[Organism] AND (23S ribosomal RNA[Title] OR rrl[Gene Name]) \
                    AND 1000:3000[SLEN]",
            default_output: "myco_rrl.fasta",
        }),
        other => anyhow::bail!(
            "Unknown gene target '{}'. Valid options: 16s, hsp65, rpob, erm41, rrl",
            other
        ),
    }
}

// ── NCBI Entrez base URL ─────────────────────────────────────────────────────

const BASE: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let args = Args::parse();
    let target = resolve_target(&args.gene)?;
    let output = args
        .output
        .clone()
        .unwrap_or_else(|| target.default_output.to_string());

    let client = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(60))
        .user_agent(format!("myco_downloader/0.1 ({})", args.email))
        .build()?;

    println!("=== Mycobacterium {} Sequence Downloader ===", target.label);
    println!("Query  : {}", target.query);
    println!("Output : {}", output);
    println!();

    // ── Step 1: esearch — get total count ────────────────────────────────────
    println!("[1/3] Searching NCBI nucleotide database...");

    let search_url = build_esearch_url(&args, target.query, 0, 1);
    let search_resp: serde_json::Value = client
        .get(&search_url)
        .send()
        .context("esearch request failed")?
        .json()
        .context("esearch JSON parse failed")?;

    let total: usize = search_resp["esearchresult"]["count"]
        .as_str()
        .unwrap_or("0")
        .parse()
        .unwrap_or(0);

    println!("  Found {:>6} sequences on NCBI", total);

    if total == 0 {
        println!("No sequences found — check query or network.");
        return Ok(());
    }

    if args.count_only {
        println!("\nCount-only mode — exiting without download.");
        return Ok(());
    }

    let fetch_total = if args.max > 0 {
        args.max.min(total)
    } else {
        total
    };
    println!("  Will download: {} sequences", fetch_total);
    println!();

    // ── Step 2: esearch with usehistory=y ────────────────────────────────────
    println!("[2/3] Storing results on NCBI history server...");
    let history_url = build_esearch_url_history(&args, target.query, fetch_total);
    let history_resp: serde_json::Value = client
        .get(&history_url)
        .send()
        .context("esearch history request failed")?
        .json()
        .context("esearch history JSON parse failed")?;

    let web_env = history_resp["esearchresult"]["webenv"]
        .as_str()
        .unwrap_or("")
        .to_string();
    let query_key = history_resp["esearchresult"]["querykey"]
        .as_str()
        .unwrap_or("1")
        .to_string();

    println!(
        "  History server ready. WebEnv: {}...",
        &web_env[..20.min(web_env.len())]
    );
    println!();

    // ── Step 3: efetch in batches ─────────────────────────────────────────────
    println!(
        "[3/3] Downloading sequences in batches of {}...",
        args.batch
    );

    let mut out = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&output)
        .context("Cannot create output file")?;

    let pb = ProgressBar::new(fetch_total as u64);
    pb.set_style(
        ProgressStyle::with_template(
            "[{elapsed_precise}] {bar:40.cyan/blue} {pos:>6}/{len:6} seqs  {msg}",
        )
        .unwrap()
        .progress_chars("=>-"),
    );

    let mut downloaded = 0usize;
    let mut retstart = 0usize;

    while retstart < fetch_total {
        let batch_size = args.batch.min(fetch_total - retstart);

        let fetch_url = build_efetch_url(&args, &web_env, &query_key, retstart, batch_size);

        let fasta_chunk = retry_get(&client, &fetch_url, 3)?;
        let seq_count = fasta_chunk.matches('>').count();

        out.write_all(fasta_chunk.as_bytes())
            .context("Write to output file failed")?;

        downloaded += seq_count;
        retstart += batch_size;

        pb.set_position(downloaded as u64);
        pb.set_message(format!("batch +{}", seq_count));

        let delay = if args.api_key.is_some() { 120 } else { 340 };
        thread::sleep(Duration::from_millis(delay));

        if seq_count == 0 {
            break;
        }
    }

    pb.finish_with_message("done");
    println!();
    println!("=== Complete ===");
    println!("  Sequences downloaded : {}", downloaded);
    println!("  Output file          : {}", output);
    println!();
    println!("Next steps:");
    println!(
        "  makeblastdb -in {} -dbtype nucl -out myco_{}_db",
        output, args.gene
    );
    println!(
        "  blastn -query query.fasta -db myco_{}_db -outfmt 6 -perc_identity 99",
        args.gene
    );

    Ok(())
}

// ── URL builders ─────────────────────────────────────────────────────────────

fn api_key_param(args: &Args) -> String {
    match &args.api_key {
        Some(k) => format!("&api_key={}", k),
        None => String::new(),
    }
}

fn build_esearch_url(args: &Args, query: &str, retstart: usize, retmax: usize) -> String {
    format!(
        "{}/esearch.fcgi?db=nucleotide&term={}&retstart={}&retmax={}\
         &retmode=json&email={}{}",
        BASE,
        urlencode(query),
        retstart,
        retmax,
        args.email,
        api_key_param(args)
    )
}

fn build_esearch_url_history(args: &Args, query: &str, retmax: usize) -> String {
    format!(
        "{}/esearch.fcgi?db=nucleotide&term={}&usehistory=y&retmax={}\
         &retmode=json&email={}{}",
        BASE,
        urlencode(query),
        retmax,
        args.email,
        api_key_param(args)
    )
}

fn build_efetch_url(
    args: &Args,
    web_env: &str,
    query_key: &str,
    retstart: usize,
    retmax: usize,
) -> String {
    format!(
        "{}/efetch.fcgi?db=nucleotide&WebEnv={}&query_key={}\
         &retstart={}&retmax={}&rettype=fasta&retmode=text&email={}{}",
        BASE,
        web_env,
        query_key,
        retstart,
        retmax,
        args.email,
        api_key_param(args)
    )
}

// ── Helpers ───────────────────────────────────────────────────────────────────

fn urlencode(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            ' ' => "%20".to_string(),
            '[' => "%5B".to_string(),
            ']' => "%5D".to_string(),
            '/' => "%2F".to_string(),
            ':' => "%3A".to_string(),
            '"' => "%22".to_string(),
            '(' => "%28".to_string(),
            ')' => "%29".to_string(),
            c => c.to_string(),
        })
        .collect()
}

fn retry_get(client: &reqwest::blocking::Client, url: &str, attempts: u32) -> Result<String> {
    for attempt in 1..=attempts {
        match client.get(url).send() {
            Ok(resp) if resp.status().is_success() => {
                return resp.text().context("Failed to read response body");
            }
            Ok(resp) => {
                eprintln!(
                    "  HTTP {} on attempt {}/{}",
                    resp.status(),
                    attempt,
                    attempts
                );
            }
            Err(e) => {
                eprintln!("  Request error on attempt {}/{}: {}", attempt, attempts, e);
            }
        }
        if attempt < attempts {
            thread::sleep(Duration::from_secs(5 * attempt as u64));
        }
    }
    anyhow::bail!("All {} attempts failed for URL: {}", attempts, url)
}
