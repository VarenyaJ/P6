"""
End-to-end unit test for DefaultMapper.apply_mapping and construct_phenopacket_for_patient.

We build tiny genotype/phenotype DataFrames for a single patient and assert:
- stats counters are set,
- a Phenopacket is produced with expected minimal content.
"""

import pandas as pd
import hpotk
from stairval.notepad import create_notepad
from P6.mapper import DefaultMapper

HPO_PATH = "tests/data/hp.v2024-04-26.json.gz"


def test_apply_mapping_builds_phenopackets_and_stats():
    mapper = DefaultMapper(hpotk.load_minimal_ontology(HPO_PATH))
    note = create_notepad("apply")

    # DataFrames use the index as the patient id; mapper will bring it into a column.
    geno = pd.DataFrame(
        {
            "contact_email": ["user@example.com"],
            "phasing": [1],
            "chromosome": ["chr16"],
            "start_position": [100],
            "end_position": [100],
            "reference": ["A"],
            "alternate": ["G"],
            "gene_symbol": ["GENE1"],
            "hgvsg": ["chr16:g.100A>G"],
            "hgvsc": ["NM_000000.0:c.100A>G"],
            "hgvsp": ["NP_000000.0:p.(Lys34Glu)"],
            "zygosity": ["het"],
            "inheritance": ["inherited"],
        },
        index=pd.Index(["P100"], name="patient"),
    )

    pheno = pd.DataFrame(
        {
            "hpo_id": ["HP:510"],  # will normalize to HP:0000510
            "date_of_observation": [20200101],
            "status": [1],
        },
        index=pd.Index(["P100"], name="patient"),
    )

    tables = {"genotype": geno, "phenotype": pheno}
    packets = mapper.apply_mapping(tables, note)

    print("ERRS", list(note.errors()))
    print("WARNS", list(note.warnings()))

    # Stats should reflect 1 genotype + 1 phenotype for 1 patient
    assert mapper.stats["genotypes"] == 1
    assert mapper.stats["phenotypes"] == 1
    assert mapper.stats["patients"] == 1
    assert len(packets) == 1

    pkt = packets[0]
    assert pkt.id == "P100"
    assert pkt.subject.id == "P100"
    # one phenotypic feature with normalized ID
    assert pkt.phenotypic_features[0].type.id == "HP:0000510"
    # one interpretation with an HGVS expression
    assert (
        pkt.interpretations[0]
        .diagnosis.genomic_interpretations[0]
        .variant_interpretation.variation_descriptor.expressions[0]
        .value
        == "16:g.100A>G"
    )

    # No mapping errors expected for this minimal happy path
    assert not note.has_errors(include_subsections=True)
