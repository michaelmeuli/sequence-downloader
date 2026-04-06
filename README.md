# myco_16s_downloader

Rust tool to download **Mycobacterium 16S rRNA** sequences from NCBI Entrez,
ready to use as a local BLAST database.

## Requirements

- Rust 1.70+ (`rustup` from https://rustup.rs)
- Internet access to NCBI

## Build

```bash
cargo build --release
```

## Usage

```bash
# Download all sequences (may be several thousand)
cargo run --release -- --output myco_16s.fasta

# Download up to 500 sequences
cargo run --release -- --output myco_16s.fasta --max 500

# With your email (required by NCBI policy) and API key (10x faster)
cargo run --release -- \
  --output myco_16s.fasta \
  --email your@email.com \
  --api-key YOUR_NCBI_API_KEY

# Just check how many sequences exist without downloading
cargo run --release -- --count-only
```

## NCBI Search Query Used

```
Mycobacterium[Organism] AND 16S ribosomal RNA[Title]
AND 400:1800[SLEN] AND biomol_rrna[PROP]
```

- `Mycobacterium[Organism]` — genus-level filter
- `16S ribosomal RNA[Title]` — gene name in record title
- `400:1800[SLEN]` — typical 16S amplicon length range
- `biomol_rrna[PROP]` — ribosomal RNA molecule type

Remove `biomol_rrna[PROP]` to also include genomic sequences containing 16S.

## Build local BLAST database

After download, run (requires NCBI BLAST+ installed):

```bash
makeblastdb -in myco_16s.fasta -dbtype nucl -out myco_16s_db \
  -title "Mycobacterium 16S rRNA"

blastn -query query.fasta -db myco_16s_db \
  -outfmt 6 -perc_identity 99 -max_target_seqs 10
```

## Output format

Standard FASTA, e.g.:
```
>NR_044908.1 Mycobacterium tuberculosis H37Rv 16S ribosomal RNA, complete sequence
AGAGTTTGATCCTGGCTCAGGACGAACGCTGGCGGCGTGCTTAACACATGCAAGTCGAACGGAAAGGTCT...
```

## CLI options

| Flag | Default | Description |
|------|---------|-------------|
| `--output` / `-o` | `myco_16s.fasta` | Output FASTA file |
| `--max` / `-m` | `0` (all) | Max sequences to download |
| `--email` / `-e` | `researcher@example.com` | Your email (NCBI policy) |
| `--api-key` / `-a` | none | NCBI API key (optional) |
| `--batch` / `-b` | `200` | Sequences per efetch call |
| `--count-only` | false | Print count only, no download |

## Notes

- Without an API key NCBI allows 3 requests/second; the tool respects this automatically
- Get a free NCBI API key at: https://www.ncbi.nlm.nih.gov/account/
- For hsp65 or rpoB, change `SEARCH_QUERY` in `src/main.rs` accordingly
