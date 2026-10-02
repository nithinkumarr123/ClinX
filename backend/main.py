from fastapi import FastAPI, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from loguru import logger
import sys

from config import config
from models.schemas import (
    VariantInput, AnalysisOutput, MutationType,
    DNAScore, ProteinScore, ProteinAnalysis,
    TranscriptInfo, UniProtInfo, ACMGResult, ClinVarInfo,
    ACMGClassification,
)
from services.ml_client import ml_client
from services.refseq_service import refseq_service
from services.transcript_service import transcript_service
from services.clinvar_service import clinvar_service
from services.calibration_service import calibration_service
from services.acmg_service import acmg_service
from services.uniprot_service import uniprot_service
from data.transcript_metadata import TRANSCRIPT_METADATA

logger.remove()
logger.add(sys.stdout, level=config.LOG_LEVEL)

app = FastAPI(title=config.API_TITLE, version=config.API_VERSION)
app.add_middleware(
    CORSMiddleware,
    allow_origins=config.CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}


def revcomp(seq: str) -> str:
    return "".join(COMPLEMENT.get(b, "N") for b in reversed(seq.upper()))


# Cache CDS/protein per gene
_cds_cache = {}


async def get_gene_sequences(gene: str) -> dict:
    if gene in _cds_cache:
        return _cds_cache[gene]

    meta = TRANSCRIPT_METADATA.get(gene)
    if not meta:
        raise HTTPException(404, f"No metadata for {gene}")

    result = await refseq_service.get_cds_and_protein(
        gene, meta["refseq_mrna"], meta["refseq_protein"]
    )
    if not result:
        raise HTTPException(502, f"Could not fetch sequences for {gene}")

    _cds_cache[gene] = result
    return result


def get_codon_at_cds_position(cds: str, cds_pos: int) -> str:
    idx = cds_pos - 1
    start = (idx // 3) * 3
    return cds[start:start + 3]


def apply_mutation(seq: str, pos: int, ref: str, alt: str) -> str:
    """Apply a substitution to `seq` at 0-based `pos`. Returns mutated seq."""
    if pos + len(ref) > len(seq):
        return seq
    actual = seq[pos:pos + len(ref)]
    if actual != ref:
        logger.error(f"❌ REF mismatch: expected '{ref}', got '{actual}'")
        return seq
    return seq[:pos] + alt + seq[pos + len(ref):]


def normalize_alleles_to_transcript(gene: str, ref: str, alt: str) -> tuple:
    """
    Convert alleles from genomic plus-strand to transcript-strand.
    For minus-strand genes, this is a complement (single base) or revcomp (multi-base).
    """
    meta = TRANSCRIPT_METADATA.get(gene)
    if not meta:
        return ref.upper(), alt.upper()

    if meta["strand"] == -1:
        if len(ref) == 1 and len(alt) == 1:
            return COMPLEMENT[ref.upper()], COMPLEMENT[alt.upper()]
        # For indels, revcomp the whole allele string
        return revcomp(ref), revcomp(alt)
    return ref.upper(), alt.upper()


def protein_fragment(protein: str, aa_pos: int, alt_aa: str = None) -> str:
    if not protein or not aa_pos or aa_pos <= 0:
        return protein[:30] + "..." if protein else ""
    aa_pos = min(aa_pos, len(protein))
    start = max(0, aa_pos - 10)
    end = min(len(protein), aa_pos + 10)
    frag = list(protein[start:end])
    idx = aa_pos - start - 1
    if 0 <= idx < len(frag):
        marker = alt_aa if alt_aa and alt_aa != frag[idx] and alt_aa != '?' else frag[idx]
        frag[idx] = f"[{marker}]"
    return ''.join(frag)


def is_frameshift(ref: str, alt: str) -> bool:
    return len(ref) != len(alt)


@app.post("/analyze", response_model=AnalysisOutput)
async def analyze_variant(variant: VariantInput):
    logger.info(f"🔬 {variant.gene} {variant.position} {variant.reference_allele}>{variant.alternate_allele}")

    # ===== 1. Fetch CDS and protein from NCBI RefSeq =====
    seq_data = await get_gene_sequences(variant.gene)
    ref_cds = seq_data["cds"]
    ref_protein = seq_data["protein"]
    transcript_id = seq_data["transcript_id"]

    logger.info(f"📊 CDS={len(ref_cds)}bp | Protein={len(ref_protein) if ref_protein else 0}aa")

    # ===== 2. Map variant to CDS/AA =====
    mapped = await transcript_service.map_variant(
        variant.gene, variant.position,
        variant.reference_allele, variant.alternate_allele,
    )
    if not mapped:
        raise HTTPException(400, f"Could not map {variant.gene}:{variant.position}")

    if mapped.get("status") != "coding":
        logger.warning(f"⚠️ Variant is {mapped['status']}")
        return AnalysisOutput(
            variant=variant.model_dump(),
            transcript=TranscriptInfo(
                transcript_id=transcript_id,
                gene_symbol=variant.gene,
                cds_sequence=ref_cds[:200],
                protein_sequence=ref_protein[:100] if ref_protein else "",
                variant_in_cds=False,
                cds_position=None,
                aa_position=None,
                source=mapped.get("source"),
                status=mapped.get("status"),
            ),
            dna_score=DNAScore(interpretation="non-coding"),
            protein_analysis=ProteinAnalysis(
                mutation_type=MutationType.INTRONIC if mapped.get("status") == "intronic" else MutationType.UTR,
                ref_aa='?', alt_aa='?', aa_position=0,
                protein_score=ProteinScore(),
                ref_fragment="", alt_fragment="",
            ),
            uniprot=UniProtInfo(protein_name=f"{variant.gene} protein", function="N/A"),
            clinvar=ClinVarInfo(significance="Non-coding variant"),
            acmg=ACMGResult(
                classification=ACMGClassification.VUS,
                score=0.0,
                evidence={"Note": [f"Variant is {mapped['status']}"]},
                summary=f"Variant is {mapped['status']} - insufficient evidence",
            ),
        )

    cds_pos = mapped["cds_position"]
    aa_pos = mapped["aa_position"]

    # ===== 3. Normalize alleles to transcript strand =====
    # The CDS from RefSeq is on the TRANSCRIPT strand.
    # Input alleles are genomic plus-strand. For minus-strand genes, flip them.
    ref_tx, alt_tx = normalize_alleles_to_transcript(
        variant.gene,
        variant.reference_allele,
        variant.alternate_allele,
    )
    logger.info(f"🧬 Alleles: input {variant.reference_allele}>{variant.alternate_allele} "
                f"→ transcript {ref_tx}>{alt_tx}")

    frameshift = is_frameshift(ref_tx, alt_tx)
    alt_cds = ref_cds
    ref_aa = '?'
    alt_aa = '?'
    mtype = MutationType.SYNONYMOUS

    # ===== 4. Mutate CDS using transcript-strand alleles =====
    if frameshift:
        mtype = MutationType.FRAMESHIFT
        i0 = cds_pos - 1
        if len(alt_tx) == 0:
            alt_cds = ref_cds[:i0] + ref_cds[i0 + len(ref_tx):]
        else:
            alt_cds = ref_cds[:i0] + alt_tx + ref_cds[i0 + len(ref_tx):]
    else:
        if 0 < cds_pos <= len(ref_cds):
            alt_cds = apply_mutation(ref_cds, cds_pos - 1, ref_tx, alt_tx)
            ref_codon = get_codon_at_cds_position(ref_cds, cds_pos)
            alt_codon = get_codon_at_cds_position(alt_cds, cds_pos)
            ref_aa = CODON_TABLE.get(ref_codon, '?')
            alt_aa = CODON_TABLE.get(alt_codon, '?')
            if ref_aa == alt_aa:
                mtype = MutationType.SYNONYMOUS
            elif alt_aa == '*':
                mtype = MutationType.NONSENSE
            else:
                mtype = MutationType.MISSENSE
            logger.info(f"✅ Codon: {ref_codon}→{alt_codon}, AA: {ref_aa}→{alt_aa} ({mtype.value})")

    # ===== 5. DNA window → Nucleotide Transformer =====
    dna_delta = 0.0
    if 0 < cds_pos <= len(ref_cds):
        anchor = cds_pos
        w = 200
        start = max(0, anchor - w // 2 - 1)
        end = min(len(ref_cds), anchor + w // 2 - 1)
        if len(alt_cds) < end:
            alt_cds = alt_cds.ljust(end, 'N')
        try:
            logger.info(f"🌐 DNA scorer: window={end-start}bp, anchor={anchor}")
            dna = await ml_client.score_dna(ref_cds[start:end], alt_cds[start:end], anchor - start - 1)
            dna_delta = dna.get("delta_score", 0.0)
            logger.info(f"✅ DNA delta: {dna_delta:.3f}")
        except Exception as e:
            logger.error(f"❌ DNA failed: {e}")

    if (frameshift or mtype == MutationType.NONSENSE) and dna_delta >= -0.5:
        dna_delta = -2.5
        logger.info(f"🔴 Forced disruptive DNA delta for {mtype.value}")

    calibrated = calibration_service.calibrate_dna_score(dna_delta, variant.gene)
    dna_score = DNAScore(
        delta_score=calibrated["calibrated_delta"],
        interpretation=calibrated["interpretation"],
    )

    # ===== 6. Protein window → ESM-2 =====
    protein_score = ProteinScore()
    if ref_protein:
        if 0 < aa_pos <= len(ref_protein):
            p_anchor = aa_pos
        elif aa_pos and aa_pos > len(ref_protein):
            p_anchor = len(ref_protein)
        else:
            p_anchor = max(1, len(ref_protein) // 2)

        if frameshift:
            ref_aa_esm = ref_protein[p_anchor - 1] if 0 < p_anchor <= len(ref_protein) else 'X'
            alt_protein_local = refseq_service.translate_cds(alt_cds) if hasattr(refseq_service, "translate_cds") else ""
            alt_aa_esm = alt_protein_local[p_anchor - 1] if 0 < p_anchor <= len(alt_protein_local) else 'X'
        else:
            ref_aa_esm = ref_aa if ref_aa != '?' else 'X'
            alt_aa_esm = alt_aa if alt_aa != '?' else 'X'

        if mtype == MutationType.NONSENSE:
            alt_aa_esm = '*'
            ref_aa_esm = ref_aa if ref_aa != '?' else 'Q'

        w = 50
        p_start = max(0, p_anchor - w)
        p_end = min(len(ref_protein), p_anchor + w)
        ref_window_p = ref_protein[p_start:p_end]

        if not ref_window_p:
            ref_window_p = ref_protein[:100]
            p_anchor = 1
            p_start = 0

        if frameshift:
            alt_window_p = (alt_protein_local[p_start:min(len(alt_protein_local), p_end)]
                            if alt_protein_local else ref_window_p)
            if len(alt_window_p) < len(ref_window_p):
                alt_window_p = alt_window_p.ljust(len(ref_window_p), 'X')
        else:
            alt_window_p = list(ref_window_p)
            local_idx = p_anchor - p_start - 1
            if 0 <= local_idx < len(alt_window_p):
                alt_window_p[local_idx] = alt_aa_esm
            alt_window_p = ''.join(alt_window_p)

        try:
            logger.info(f"🌐 ESM-2 scorer: window={len(ref_window_p)}aa, {ref_aa_esm}→{alt_aa_esm}")
            esm = await ml_client.score_protein(
                ref_window_p, p_anchor - p_start, ref_aa_esm, alt_aa_esm
            )
            delta = esm.get("delta_score", 0.0)
            impact = esm.get("impact", "tolerated")
            if frameshift or mtype == MutationType.NONSENSE:
                if delta > -2.0:
                    delta = min(delta, -5.0)
                    impact = "deleterious"
            protein_score = ProteinScore(delta_score=delta, impact=impact)
            logger.info(f"✅ ESM-2 delta: {delta:.3f}, impact: {impact}")
        except Exception as e:
            logger.error(f"❌ Protein failed: {e}")
            if frameshift or mtype == MutationType.NONSENSE:
                protein_score = ProteinScore(delta_score=-5.0, impact="deleterious")

    # ===== 7. ClinVar =====
    clinvar = await clinvar_service.get_clinvar_for_variant(
        variant.gene, variant.reference_allele, variant.alternate_allele, variant.position
    )

    # ===== 8. ACMG =====
    classification, posterior, evidence, summary = acmg_service.classify(
        mtype, dna_score.delta_score, protein_score.delta_score, clinvar
    )

    uniprot_info = await uniprot_service.get_protein_info(variant.gene)

    return AnalysisOutput(
        variant=variant.model_dump(),
        transcript=TranscriptInfo(
            transcript_id=transcript_id,
            gene_symbol=variant.gene,
            cds_sequence=ref_cds[:200],
            protein_sequence=ref_protein[:100] if ref_protein else "",
            variant_in_cds=True,
            cds_position=cds_pos,
            aa_position=aa_pos,
            source=mapped.get("source"),
            status=mapped.get("status"),
        ),
        dna_score=dna_score,
        protein_analysis=ProteinAnalysis(
            mutation_type=mtype,
            ref_aa=ref_aa,
            alt_aa=alt_aa,
            aa_position=aa_pos or 0,
            protein_score=protein_score,
            ref_fragment=protein_fragment(ref_protein, aa_pos),
            alt_fragment=protein_fragment(ref_protein, aa_pos, alt_aa if alt_aa != '?' else None),
        ),
        uniprot=UniProtInfo(**uniprot_info),
        clinvar=ClinVarInfo(
            significance=clinvar.get("clinical_significance", "Unknown"),
            conditions=clinvar.get("conditions", ["Not specified"]),
            review_status=clinvar.get("review_status", "No review"),
            found=clinvar.get("found", False),
        ),
        acmg=ACMGResult(
            classification=classification,
            score=float(posterior),
            evidence=evidence,
            summary=summary,
        ),
    )


@app.get("/health")
async def health():
    ml_health = await ml_client.health_check()
    return {"status": "ok", "version": config.API_VERSION, "ml_server": ml_health}


if __name__ == "__main__":
    import uvicorn
    uvicorn.run(app, host=config.API_HOST, port=config.API_PORT)