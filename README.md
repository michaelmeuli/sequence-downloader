# myco_downloader

Rust tool to download **Mycobacteriaceae gene sequences** from NCBI Entrez,
ready to use as local BLAST databases.

## Supported targets

| Gene | Label | Notes |
|------|-------|-------|
| `rrs` | 16S rRNA | type material, 400–3000 bp |
| `hsp65` | hsp65 / groEL2 | type material, 400–3000 bp |
| `rpob` | rpoB | type material, 400–3000 bp |
| `erm41` | erm(41) | deduplication applied automatically |
| `rrl` | 23S rRNA | type material, 400–3000 bp |
| `all` | all of the above | each written to its own file |

## Requirements

- Rust 1.70+ (`rustup` from https://rustup.rs)
- Internet access to NCBI

## Build

```bash
cargo build --release
```

## Usage

```bash
# Download 16S rRNA sequences (default)
cargo run --release

# Download a specific gene
cargo run --release -- --gene hsp65

# Download all genes at once
cargo run --release -- --gene all

# Limit to 500 sequences
cargo run --release -- --gene rrs --max 500

# With your email (required by NCBI policy) and API key (10x faster)
cargo run --release -- --gene all \
  --email your@email.com \
  --api-key YOUR_NCBI_API_KEY

# Check sequence counts without downloading
cargo run --release -- --gene all --count-only
```

## Build local BLAST databases

After download, run (requires NCBI BLAST+ installed):

```bash
makeblastdb -in myco_rrs.fasta -dbtype nucl -out myco_rrs_db
blastn -query query.fasta -db myco_rrs_db -outfmt 6 -perc_identity 99
```

## Output

Standard FASTA files, one per gene target (e.g. `myco_rrs.fasta`):

```
>NR_044908.1 Mycobacterium tuberculosis H37Rv 16S ribosomal RNA, complete sequence
AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGTCT...
```

## CLI options

| Flag | Default | Description |
|------|---------|-------------|
| `--gene` / `-g` | `rrs` | Target gene (`rrs`, `hsp65`, `rpob`, `erm41`, `rrl`, `all`) |
| `--output` / `-o` | `myco_<gene>.fasta` | Output file (ignored for `--gene all`) |
| `--max` / `-m` | `0` (all) | Max sequences to download |
| `--email` / `-e` | `researcher@example.com` | Your email (NCBI policy) |
| `--api-key` / `-a` | none | NCBI API key (optional) |
| `--batch` / `-b` | `200` | Sequences per efetch call |
| `--count-only` | false | Print count only, no download |

## Notes

- Without an API key NCBI allows 3 requests/second; the tool respects this automatically
- Get a free NCBI API key at: https://www.ncbi.nlm.nih.gov/account/
- `erm41` sequences are automatically deduplicated by accession after download
