/// myco_16s_downloader
///
/// Downloads 16S rRNA sequences for Mycobacterium species from NCBI
/// using the Entrez E-utilities API (esearch + efetch).
///
/// Usage:
///   cargo run -- --output myco_16s.fasta
///   cargo run -- --output myco_16s.fasta --max 500 --email your@email.com
///
/// NCBI Entrez API endpoints used:
///   esearch: search nucleotide DB and get list of GI/UIDs
///   efetch:  fetch sequences in FASTA format by UID batch

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
    name = "myco_16s_downloader",
    about = "Download Mycobacterium 16S rRNA sequences from NCBI Entrez"
)]
struct Args {
    /// Output FASTA file path
    #[arg(short, long, default_value = "myco_16s.fasta")]
    output: String,

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

// ── NCBI Entrez base URL ─────────────────────────────────────────────────────

const BASE: &str = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils";

// ── Search query ─────────────────────────────────────────────────────────────
//
// Targets: 16S rRNA gene sequences from genus Mycobacterium
// - "Mycobacterium"[Organism]   → restrict to genus
// - 16S rRNA[Title]             → gene name in record title
// - 400:1800[SLEN]              → typical 16S amplicon length range
// - biomol_rrna[PROP]           → ribosomal RNA molecule type
// - srcdb_refseq[PROP]          → RefSeq entries only (curated)
//   (remove srcdb_refseq to get all GenBank entries — much larger set)

const SEARCH_QUERY: &str =
    "Mycobacterium[Organism] AND 16S ribosomal RNA[Title] \
     AND 400:1800[SLEN] AND biomol_rrna[PROP]";

// ── Main ─────────────────────────────────────────────────────────────────────

fn main() -> Result<()> {
    let args = Args::parse();
    let client = reqwest::blocking::Client::builder()
        .timeout(Duration::from_secs(60))
        .user_agent(format!("myco_16s_downloader/0.1 ({})", args.email))
        .build()?;

    println!("=== Mycobacterium 16S rRNA Sequence Downloader ===");
    println!("Query  : {}", SEARCH_QUERY);
    println!("Output : {}", args.output);
    println!();

    // ── Step 1: esearch — get total count and WebEnv/query_key ───────────────
    println!("[1/3] Searching NCBI nucleotide database...");

    let search_url = build_esearch_url(&args, 0, 1);
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

    let web_env = search_resp["esearchresult"]["webenv"]
        .as_str()
        .unwrap_or("")
        .to_string();
    let query_key = search_resp["esearchresult"]["querykey"]
        .as_str()
        .unwrap_or("1")
        .to_string();

    println!("  Found {:>6} sequences on NCBI", total);

    if total == 0 {
        println!("No sequences found — check query or network.");
        return Ok(());
    }

    if args.count_only {
        println!("\nCount-only mode — exiting without download.");
        return Ok(());
    }

    // Determine how many to actually fetch
    let fetch_total = if args.max > 0 {
        args.max.min(total)
    } else {
        total
    };

    println!("  Will download: {} sequences", fetch_total);
    println!();

    // ── Step 2: Use history server (WebEnv) to batch-fetch ───────────────────
    // Re-run esearch with usehistory=y to store results server-side,
    // then efetch in batches referencing WebEnv + query_key.

    println!("[2/3] Storing results on NCBI history server...");
    let history_url = build_esearch_url_history(&args, fetch_total);
    let history_resp: serde_json::Value = client
        .get(&history_url)
        .send()
        .context("esearch history request failed")?
        .json()
        .context("esearch history JSON parse failed")?;

    let web_env2 = history_resp["esearchresult"]["webenv"]
        .as_str()
        .unwrap_or(&web_env)
        .to_string();
    let query_key2 = history_resp["esearchresult"]["querykey"]
        .as_str()
        .unwrap_or(&query_key)
        .to_string();

    println!("  History server ready. WebEnv: {}...", &web_env2[..20.min(web_env2.len())]);
    println!();

    // ── Step 3: efetch in batches ─────────────────────────────────────────────
    println!("[3/3] Downloading sequences in batches of {}...", args.batch);

    // Create/truncate output file
    let mut out = OpenOptions::new()
        .write(true)
        .create(true)
        .truncate(true)
        .open(&args.output)
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

    while downloaded < fetch_total {
        let batch_size = args.batch.min(fetch_total - downloaded);

        let fetch_url = build_efetch_url(
            &args,
            &web_env2,
            &query_key2,
            retstart,
            batch_size,
        );

        // Retry loop (NCBI can 429 or 503 occasionally)
        let fasta_chunk = retry_get(&client, &fetch_url, 3)?;

        // Count sequences fetched in this batch
        let seq_count = fasta_chunk.matches('>').count();

        out.write_all(fasta_chunk.as_bytes())
            .context("Write to output file failed")?;

        downloaded += seq_count;
        retstart += batch_size;

        pb.set_position(downloaded as u64);
        pb.set_message(format!("batch +{}", seq_count));

        // NCBI rate limit: 3 req/s without API key, 10/s with
        let delay = if args.api_key.is_some() { 120 } else { 340 };
        thread::sleep(Duration::from_millis(delay));

        if seq_count == 0 {
            break; // No more sequences
        }
    }

    pb.finish_with_message("done");
    println!();
    println!("=== Complete ===");
    println!("  Sequences downloaded : {}", downloaded);
    println!("  Output file          : {}", args.output);
    println!();
    println!("Next steps:");
    println!("  makeblastdb -in {} -dbtype nucl -out myco_16s_db", args.output);
    println!("  blastn -query query.fasta -db myco_16s_db -outfmt 6 -perc_identity 99");

    Ok(())
}

// ── URL builders ─────────────────────────────────────────────────────────────

fn api_key_param(args: &Args) -> String {
    match &args.api_key {
        Some(k) => format!("&api_key={}", k),
        None => String::new(),
    }
}

/// Initial esearch to get total count
fn build_esearch_url(args: &Args, retstart: usize, retmax: usize) -> String {
    format!(
        "{}/esearch.fcgi?db=nucleotide&term={}&retstart={}&retmax={}\
         &retmode=json&email={}{}",
        BASE,
        urlencode(SEARCH_QUERY),
        retstart,
        retmax,
        args.email,
        api_key_param(args)
    )
}

/// esearch with usehistory=y — stores results server-side for batch efetch
fn build_esearch_url_history(args: &Args, retmax: usize) -> String {
    format!(
        "{}/esearch.fcgi?db=nucleotide&term={}&usehistory=y&retmax={}\
         &retmode=json&email={}{}",
        BASE,
        urlencode(SEARCH_QUERY),
        retmax,
        args.email,
        api_key_param(args)
    )
}

/// efetch using WebEnv history key — returns FASTA text
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

/// Very simple URL encoder for the query string
fn urlencode(s: &str) -> String {
    s.chars()
        .map(|c| match c {
            ' ' => "%20".to_string(),
            '[' => "%5B".to_string(),
            ']' => "%5D".to_string(),
            '/' => "%2F".to_string(),
            ':' => "%3A".to_string(),
            '"' => "%22".to_string(),
            c => c.to_string(),
        })
        .collect()
}

/// GET with simple retry on transient errors
fn retry_get(
    client: &reqwest::blocking::Client,
    url: &str,
    attempts: u32,
) -> Result<String> {
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
