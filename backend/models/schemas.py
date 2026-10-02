from pydantic import BaseModel, Field, validator
from typing import Optional, List, Dict, Any
from enum import Enum
from config import config


class GenomeBuild(str, Enum):
    HG38 = "GRCh38"
    HG19 = "GRCh37"


class MutationType(str, Enum):
    SYNONYMOUS = "synonymous"
    MISSENSE = "missense"
    NONSENSE = "nonsense"
    FRAMESHIFT = "frameshift"
    INTRONIC = "intronic"
    UTR = "utr"


class ACMGClassification(str, Enum):
    PATHOGENIC = "Pathogenic"
    LIKELY_PATHOGENIC = "Likely Pathogenic"
    VUS = "Variant of Uncertain Significance"
    LIKELY_BENIGN = "Likely Benign"
    BENIGN = "Benign"


class VariantInput(BaseModel):
    chromosome: str
    position: int
    reference_allele: str
    alternate_allele: str
    gene: str
    genome_build: GenomeBuild = GenomeBuild.HG38

    @validator('chromosome')
    def fmt_chrom(cls, v):
        return v if v.startswith('chr') else f'chr{v}'

    @validator('reference_allele', 'alternate_allele', pre=True)
    def validate_allele(cls, v):
        if v is None:
            return ""
        v = str(v).upper()
        if v and not all(c in 'ACGT' for c in v):
            raise ValueError(f'Invalid allele: {v}')
        return v

    @validator('gene')
    def normalize_gene(cls, v):
        return config.GENE_ALIASES.get(v.upper(), v.upper())


class TranscriptInfo(BaseModel):
    transcript_id: str
    gene_symbol: str
    cds_sequence: str
    protein_sequence: str
    variant_in_cds: bool
    cds_position: Optional[int] = None
    aa_position: Optional[int] = None
    source: Optional[str] = None
    status: Optional[str] = None


class DNAScore(BaseModel):
    ref_pll: float = 0.0
    alt_pll: float = 0.0
    delta_score: float = 0.0
    interpretation: str = "neutral"


class ProteinScore(BaseModel):
    ref_score: float = 0.0
    alt_score: float = 0.0
    delta_score: float = 0.0
    impact: str = "tolerated"


class ProteinAnalysis(BaseModel):
    mutation_type: MutationType
    ref_aa: str
    alt_aa: str
    aa_position: int
    protein_score: ProteinScore
    ref_fragment: str
    alt_fragment: str


class UniProtInfo(BaseModel):
    protein_name: str
    function: str
    tissue_expression: List[str] = []


class ClinVarInfo(BaseModel):
    significance: str = "Unknown"
    conditions: List[str] = []
    review_status: str = "No review"
    found: bool = False


class ACMGResult(BaseModel):
    classification: ACMGClassification
    score: float
    evidence: Dict[str, List[str]] = {}
    summary: str


class AnalysisOutput(BaseModel):
    variant: Dict[str, Any]
    transcript: TranscriptInfo
    dna_score: DNAScore
    protein_analysis: ProteinAnalysis
    uniprot: UniProtInfo
    clinvar: Optional[ClinVarInfo] = None
    acmg: ACMGResult