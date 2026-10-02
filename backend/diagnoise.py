# save as diagnose.py in the backend folder
import asyncio
from services.sequence_service import sequence_service
from services.ucsc_service import ucsc_service
from data.transcript_metadata import TRANSCRIPT_METADATA


async def diagnose():
    gene = "BRCA1"
    meta = TRANSCRIPT_METADATA[gene]
    cds_start = meta["cds_start"]     # 43124027
    cds_end = meta["cds_end"]         # 43045802
    strand = meta["strand"]           # -1
    chrom = f"chr{meta['chromosome']}"
    cds_lo = min(cds_start, cds_end)
    cds_hi = max(cds_start, cds_end)

    sorted_exons = sorted(meta["exons"], key=lambda e: e["start"], reverse=(strand == -1))

    print(f"CDS span: {cds_lo} - {cds_hi}")
    print(f"Strand: {strand}")
    print(f"Exons in transcript order (strand={strand}):")
    total = 0
    for i, exon in enumerate(sorted_exons):
        e_start = max(exon["start"], cds_lo)
        e_end = min(exon["end"], cds_hi)
        if e_start >= e_end:
            print(f"  [{i}] SKIPPED {exon['start']}-{exon['end']} (outside CDS)")
            continue
        length = e_end - e_start + 1
        total += length
        print(f"  [{i}] {e_start}-{e_end}  ({length} bp)")

    print(f"\nTotal assembled: {total} bp (expected 5592 bp)")
    print(f"Missing: {5592 - total} bp")


asyncio.run(diagnose())