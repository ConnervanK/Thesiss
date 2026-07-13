#!/usr/bin/env python3
"""
read_pdf.py — Lightweight PDF text extractor for AI assistants.

Usage examples:
  python read_pdf.py paper.pdf                     # Extract all text
  python read_pdf.py paper.pdf --max 5             # Extract only the first 5 pages
  python read_pdf.py paper.pdf --start 3 --max 2   # Extract pages 3 and 4
  python read_pdf.py paper.pdf --search "method"   # Only print pages containing "method"
"""

import argparse
import sys
from pathlib import Path
from pypdf import PdfReader
from pypdf.errors import PdfReadError


def extract_pdf_text(
    pdf_path: str,
    start_page: int = 1,
    max_pages: int | None = None,
    search_term: str | None = None,
) -> None:
    path = Path(pdf_path).resolve()

    if not path.exists():
        print(f"[Error] File not found: {path}", file=sys.stderr)
        sys.exit(1)

    try:
        reader = PdfReader(path)
    except PdfReadError as e:
        print(f"[Error] Failed to read PDF: {e}", file=sys.stderr)
        sys.exit(1)

    total_pages = len(reader.pages)
    start_idx = max(0, start_page - 1)
    end_idx = (
        total_pages
        if max_pages is None
        else min(total_pages, start_idx + max_pages)
    )

    print(
        f"--- Document: {path.name} (Pages {start_idx + 1} to {end_idx} of {total_pages}) ---"
    )

    pages_output = 0
    for idx in range(start_idx, end_idx):
        page = reader.pages[idx]
        text = page.extract_text() or ""

        # Optional filtering to save LLM tokens during targeted queries
        if search_term and search_term.lower() not in text.lower():
            continue

        print(f"\n[--- Page {idx + 1} ---]\n")
        print(text.strip())
        pages_output += 1

    if search_term and pages_output == 0:
        print(
            f"\n[No pages contained the search term: '{search_term}']",
            file=sys.stderr,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract text from a PDF for LLM workflows."
    )
    parser.add_argument("pdf_path", help="Path to the PDF file.")
    parser.add_argument(
        "--start",
        type=int,
        default=1,
        help="1-indexed starting page (default: 1).",
    )
    parser.add_argument(
        "--max",
        type=int,
        default=None,
        help="Maximum number of pages to extract.",
    )
    parser.add_argument(
        "--search", type=str, default=None, help="Only output pages containing this term."
    )

    args = parser.parse_args()
    extract_pdf_text(args.pdf_path, args.start, args.max, args.search)