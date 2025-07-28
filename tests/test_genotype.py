import pytest
from P6.genotype import Genotype


def test_valid_genotype_instantiation():
    """A valid Genotype should be created without error."""
    g = Genotype(
        genotype_patient_ID="PAT123",
        contact_email="alice.bob+test@example.com",
        phasing=True,
        chromosome="hgvs",
        start_position=100,
        end_position=200,
        reference="A",
        alternate="T",
        gene_symbol="BRCA1",
        hgvsg="g.100A>T",
        hgvsc="c.200A>T",
        hgvsp="p.Lys67Asn",
        zygosity="heterozygous",
        inheritance="inherited"
    )
    assert isinstance(g, Genotype)


@pytest.mark.parametrize("bad_id", ["", "123-!"])
def test_invalid_patient_id_raises(bad_id):
    """Non‑alphanumeric IDs must trigger a ValueError."""
    with pytest.raises(ValueError):
        Genotype(
            genotype_patient_ID=bad_id,
            contact_email="foo@bar.com",
            phasing=False,
            chromosome="ucsc",
            start_position=0,
            end_position=1,
            reference="C",
            alternate="G",
            gene_symbol="GENE",
            hgvsg="g.1C>G",
            hgvsc="c.1C>G",
            hgvsp="p.Ala1Gly",
            zygosity="homozygous",
            inheritance="unknown"
        )


@pytest.mark.parametrize("bad_email", ["noatsymbol", "foo@bar", "@foo.com"])
def test_invalid_email_format_raises(bad_email):
    """Malformed emails must trigger a ValueError."""
    with pytest.raises(ValueError):
        Genotype(
            genotype_patient_ID="PAT1",
            contact_email=bad_email,
            phasing=False,
            chromosome="refseq",
            start_position=10,
            end_position=20,
            reference="G",
            alternate="A",
            gene_symbol="GENE2",
            hgvsg="g.10G>A",
            hgvsc="c.10G>A",
            hgvsp="p.Gly4Asp",
            zygosity="mosaic",
            inheritance="de_novo_mutation"
        )
