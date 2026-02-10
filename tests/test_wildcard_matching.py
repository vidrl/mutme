import textwrap

from mutme.mutation import compare_nextclade_to_annotations


def _write(p, content: str) -> None:
    p.write_text(textwrap.dedent(content).lstrip("\n"), encoding="utf-8")


def _minimal_nextclade_tsv(
    *, subs: str = "", ins: str = "", dels: str = "", stops: str = "", qc: str = "good"
) -> str:
    """
    Build a single-row Nextclade TSV with required headers.
    Values must be comma-separated mutation tokens (or empty).
    """
    return f"""\
    seqName\tqc.overallStatus\taaSubstitutions\taaDeletions\taaInsertions\tqc.stopCodons.stopCodons
    sample_01\t{qc}\t{subs}\t{dels}\t{ins}\t{stops}
    """


def test_substitution_wildcard_matches_and_merges_with_exact(tmp_path):
    ann = tmp_path / "ann.csv"
    nxt = tmp_path / "next.tsv"

    # Both wildcard (S:N87X) and exact (S:N87Y) should match detected S:N87Y,
    # and Option B should emit ONE hit row with merged fields.
    _write(
        ann,
        """
        mutation,phenotype,drug,comment
        S:N87X,resistant,aciclovir,from wildcard
        S:N87Y,,valganciclovir,from exact
        """,
    )
    _write(
        nxt,
        _minimal_nextclade_tsv(subs="S:N87Y"),
    )

    reports = compare_nextclade_to_annotations(
        nxt,
        ann,
        allow_x_wildcards=True,
    )

    assert len(reports) == 1
    r = reports[0]
    assert r.seq_name == "sample_01"
    assert r.qc_status == "good"

    assert len(r.hits) == 1
    hit = r.hits[0]

    # mutation column should be the detected mutation
    assert hit.mutation == "S:N87Y"

    # phenotype merges from wildcard, drug merges from wildcard+exact
    assert hit.annotations["phenotype"] == "resistant"
    assert hit.annotations["drug"] in [
        "aciclovir; valganciclovir",
        "valganciclovir; aciclovir",
    ]  # NOTE MERGE LOGIC NOT DETERMINISTIC

    # comments merged (order based on merge order; we just check both appear)
    assert "from wildcard" in hit.comments
    assert "from exact" in hit.comments


def test_insertion_wildcard_matches_and_merges_with_exact(tmp_path):
    ann = tmp_path / "ann.csv"
    nxt = tmp_path / "next.tsv"

    _write(
        ann,
        """
        mutation,tag,comment
        S:345:NXY,marker,from ins wildcard
        S:345:NQY,specific,from exact ins
        """,
    )
    _write(
        nxt,
        _minimal_nextclade_tsv(ins="S:345:NQY"),
    )

    reports = compare_nextclade_to_annotations(
        nxt,
        ann,
        allow_x_wildcards=True,
    )

    assert len(reports) == 1
    hit = reports[0].hits[0]
    assert hit.mutation == "S:345:NQY"

    # merged tag field
    assert hit.annotations["tag"] in [
        "marker; specific",
        "specific; marker",
    ]  # NOTE MERGE LOGIC NOT DETERMINISTIC
    assert "from ins wildcard" in hit.comments
    assert "from exact ins" in hit.comments


def test_wildcards_disabled_preserves_exact_only(tmp_path):
    ann = tmp_path / "ann.csv"
    nxt = tmp_path / "next.tsv"

    _write(
        ann,
        """
        mutation,phenotype
        S:N87X,resistant
        """,
    )
    _write(
        nxt,
        _minimal_nextclade_tsv(subs="S:N87Y"),
    )

    reports = compare_nextclade_to_annotations(
        nxt,
        ann,
        allow_x_wildcards=False,
    )

    assert len(reports) == 1
    assert reports[0].hits == ()  # no exact match => no hits


def test_insertion_wildcard_is_length_sensitive(tmp_path):
    ann = tmp_path / "ann.csv"
    nxt = tmp_path / "next.tsv"

    # Annotation expects 3-AA insertion string with middle wildcard: N X Y
    _write(
        ann,
        """
        mutation,tag
        S:345:NXY,marker
        """,
    )

    # Detected insertion is different length => should NOT match
    _write(
        nxt,
        _minimal_nextclade_tsv(ins="S:345:NQYY"),
    )

    reports = compare_nextclade_to_annotations(
        nxt,
        ann,
        allow_x_wildcards=True,
    )

    assert len(reports) == 1
    assert reports[0].hits == ()


def test_substitution_wildcard_charset_controls_matching(tmp_path):
    ann = tmp_path / "ann.csv"
    nxt = tmp_path / "next.tsv"

    # Wildcard should match only characters in x_charset.
    _write(
        ann,
        """
        mutation,flag
        S:N87X,yes
        """,
    )
    _write(
        nxt,
        _minimal_nextclade_tsv(subs="S:N87Z"),
    )

    # Default charset is canonical 20 AAs; Z not included => no match
    reports = compare_nextclade_to_annotations(
        nxt,
        ann,
        allow_x_wildcards=True,
    )
    assert reports[0].hits == ()

    # If we allow Z explicitly, it should match
    reports2 = compare_nextclade_to_annotations(
        nxt,
        ann,
        allow_x_wildcards=True,
        x_charset="ACDEFGHIKLMNPQRSTVWYZ",
    )
    assert len(reports2[0].hits) == 1
    assert reports2[0].hits[0].mutation == "S:N87Z"
    assert reports2[0].hits[0].annotations["flag"] == "yes"


def test_multiple_wildcards_in_insertion(tmp_path):
    ann = tmp_path / "ann.csv"
    nxt = tmp_path / "next.tsv"

    _write(
        ann,
        """
        mutation,tag
        S:10:XX,any_two
        """,
    )
    _write(
        nxt,
        _minimal_nextclade_tsv(ins="S:10:QA"),
    )

    reports = compare_nextclade_to_annotations(
        nxt,
        ann,
        allow_x_wildcards=True,
    )

    assert len(reports) == 1
    assert len(reports[0].hits) == 1
    hit = reports[0].hits[0]
    assert hit.mutation == "S:10:QA"
    assert hit.annotations["tag"] == "any_two"
