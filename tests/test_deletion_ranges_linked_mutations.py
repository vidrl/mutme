from __future__ import annotations

from pathlib import Path

import pytest

from mutme.mutation import (
    compare_nextclade_to_annotations,
    get_annotation_fields,
    load_annotation_table,
)


def _write_text(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def _make_nextclade_tsv(
    path: Path,
    *,
    rows: list[dict[str, str]],
) -> Path:
    """
    Minimal Nextclade TSV writer with required columns for compare_nextclade_to_annotations.

    Notes:
    - All mutation-list columns are comma-separated within the cell, consistent with Nextclade.
    """
    header = [
        "seqName",
        "qc.overallStatus",
        "qc.stopCodons.stopCodons",
        "aaSubstitutions",
        "aaDeletions",
        "aaInsertions",
    ]
    lines = ["\t".join(header)]
    for r in rows:
        lines.append("\t".join(r.get(h, "") for h in header))
    return _write_text(path, "\n".join(lines) + "\n")


# -----------------------------------------------------------------------------
# Existing baseline tests (updated for linked_mutations separator = ';')
# -----------------------------------------------------------------------------


def test_load_annotation_table_expands_deletion_ranges_and_autolinks_by_default(
    tmp_path: Path,
) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,comment\nS:del87-89,range_hit,needs_all\n",
    )

    d = load_annotation_table(anno, delimiter=",")  # default require_full_deletion_ranges=True

    # Expansion creates canonical annotation deletion keys
    assert "S:87-" in d
    assert "S:88-" in d
    assert "S:89-" in d
    assert "S:del87-89" not in d  # range key is a convenience, not retained

    # Same annotations copied onto each expanded key
    assert d["S:87-"]["phenotype"] == "range_hit"
    assert d["S:88-"]["phenotype"] == "range_hit"
    assert d["S:89-"]["phenotype"] == "range_hit"
    assert d["S:87-"]["comment"] == "needs_all"

    # Auto-linked: each expanded deletion requires the full set (all-or-nothing)
    linked_87 = set(d["S:87-"].get("linked_mutations", "").split(";")) - {""}
    assert linked_87 == {"S:87-", "S:88-", "S:89-"}


def test_load_annotation_table_expands_deletion_ranges_without_autolink_when_disabled(
    tmp_path: Path,
) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del87-89,range_hit\n",
    )

    d = load_annotation_table(
        anno,
        delimiter=",",
        require_full_deletion_ranges=False,
    )

    # Still expands
    assert "S:87-" in d and "S:88-" in d and "S:89-" in d

    # But no implicit linked_mutations
    assert d["S:87-"].get("linked_mutations", "").strip() == ""
    assert d["S:88-"].get("linked_mutations", "").strip() == ""
    assert d["S:89-"].get("linked_mutations", "").strip() == ""


def test_compare_suppresses_partial_range_hit_by_default(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del87-89,range_hit\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "only_one_del",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "",
                "aaDeletions": "S:N87-",
                "aaInsertions": "",
            }
        ],
    )

    reports = compare_nextclade_to_annotations(
        nxt,
        anno,
        delimiter=",",
        # default require_full_deletion_ranges=True
    )
    assert len(reports) == 1
    assert reports[0].seq_name == "only_one_del"
    assert reports[0].hits == tuple()  # suppressed by auto linked_mutations gating


def test_compare_allows_partial_range_hit_when_full_range_requirement_disabled(
    tmp_path: Path,
) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del87-89,range_hit\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "only_one_del",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "",
                "aaDeletions": "S:N87-",
                "aaInsertions": "",
            }
        ],
    )

    reports = compare_nextclade_to_annotations(
        nxt,
        anno,
        delimiter=",",
        require_full_deletion_ranges=False,
    )
    assert len(reports) == 1
    r = reports[0]
    muts = {h.mutation for h in r.hits}
    assert muts == {"S:87-"}  # now allowed (no auto-linking)


def test_compare_emits_range_hits_when_all_deletions_present(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del87-89,range_hit\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "all_three_del",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "",
                "aaDeletions": "S:N87-,S:N88-,S:N89-",
                "aaInsertions": "",
            }
        ],
    )

    # Works in default mode (full-range required)...
    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert len(reports) == 1
    r = reports[0]
    muts = {h.mutation for h in r.hits}
    assert muts == {"S:87-", "S:88-", "S:89-"}

    # ...and also when disabled (it should not reduce matches)
    reports2 = compare_nextclade_to_annotations(
        nxt,
        anno,
        delimiter=",",
        require_full_deletion_ranges=False,
    )
    muts2 = {h.mutation for h in reports2[0].hits}
    assert muts2 == {"S:87-", "S:88-", "S:89-"}

    # linked_mutations is a control field, so it must not appear in annotations output
    for h in r.hits:
        assert "linked_mutations" not in h.annotations
        assert h.annotations.get("phenotype") == "range_hit"


def test_linked_mutations_gating_for_substitutions(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,S:K417N\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "missing_link",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K",
                "aaDeletions": "",
                "aaInsertions": "",
            },
            {
                "seqName": "has_link",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N",
                "aaDeletions": "",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert [r.seq_name for r in reports] == ["missing_link", "has_link"]

    assert reports[0].hits == tuple()  # gated off

    assert len(reports[1].hits) == 1
    hit = reports[1].hits[0]
    assert hit.mutation == "S:E484K"
    assert hit.annotations.get("phenotype") == "escape"
    assert "linked_mutations" not in hit.annotations


def test_linked_mutations_union_on_duplicate_rows(tmp_path: Path) -> None:
    """
    If a mutation appears multiple times with different linked_mutations, loader unions them.
    That means the hit requires all unioned linked tokens to be present.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,S:K417N\nS:E484K,escape,S:N501Y\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "has_one_link_only",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N",
                "aaDeletions": "",
                "aaInsertions": "",
            },
            {
                "seqName": "has_both_links",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N,S:N501Y",
                "aaDeletions": "",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()  # missing N501Y due to unioned requirement

    assert len(reports[1].hits) == 1
    assert reports[1].hits[0].mutation == "S:E484K"


def test_get_annotation_fields_excludes_control_columns(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,comment,linked_mutations,drug\nS:E484K,escape,foo,S:K417N,bar\n",
    )
    fields = get_annotation_fields(anno, delimiter=",")
    assert fields == ["phenotype", "drug"]  # preserves file order, excludes control cols


def test_load_annotation_table_rejects_invalid_deletion_range(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del89-87,bad\n",
    )
    with pytest.raises(Exception):  # avoid coupling on exact exception import path
        load_annotation_table(anno, delimiter=",")


# -----------------------------------------------------------------------------
# Gotcha tests: linked_mutations parsing and behavior
# -----------------------------------------------------------------------------


def test_linked_mutations_trims_whitespace_and_still_gates(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,  S:K417N  ;   S:N501Y  \n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "missing_one",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N",
                "aaDeletions": "",
                "aaInsertions": "",
            },
            {
                "seqName": "has_both",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N,S:N501Y",
                "aaDeletions": "",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()
    assert len(reports[1].hits) == 1
    assert reports[1].hits[0].mutation == "S:E484K"


def test_linked_mutations_ignores_empty_tokens_and_trailing_separators(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,;S:K417N;;S:N501Y;\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "has_both",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N,S:N501Y",
                "aaDeletions": "",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert len(reports[0].hits) == 1
    assert reports[0].hits[0].mutation == "S:E484K"


def test_linked_mutations_duplicate_tokens_do_not_break_gating(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,S:K417N;S:K417N\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "has_link",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N",
                "aaDeletions": "",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert len(reports[0].hits) == 1
    assert reports[0].hits[0].mutation == "S:E484K"


def test_linked_mutations_case_sensitive_gotcha(tmp_path: Path) -> None:
    """
    Current behavior is strict string matching. Lowercased linked tokens do not match.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,s:k417n\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "has_real_mut",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N",
                "aaDeletions": "",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()  # does not match due to case sensitivity


def test_linked_mutations_deletion_format_mismatch_gotcha(tmp_path: Path) -> None:
    """
    Deletions in the detected universe are canonicalized to {gene}:{pos}- (e.g. S:88-).
    If linked_mutations uses Nextclade-style {gene}:{aa}{pos}- (e.g. S:N88-), gating fails.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\n"
        "S:87-,del_hit,S:N88-\n",  # GOTCHA: linked token is not canonical
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "has_both_dels_raw",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "",
                "aaDeletions": "S:N87-,S:N88-",
                "aaInsertions": "",
            }
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()  # fails gating


# -----------------------------------------------------------------------------
# Gotcha tests: duplicate-row merge details
# -----------------------------------------------------------------------------


def test_linked_mutations_union_on_duplicate_rows_has_stable_serialization(tmp_path: Path) -> None:
    """
    Ensure loader stores a deterministic, de-duplicated linked_mutations string using ';'.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,S:N501Y\nS:E484K,escape,S:K417N\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    lm = d["S:E484K"].get("linked_mutations", "")
    assert lm == "S:K417N;S:N501Y"  # stable sorted serialization


def test_duplicate_rows_one_empty_linked_mutations_does_not_remove_requirement(
    tmp_path: Path,
) -> None:
    """
    Empty linked_mutations cells are ignored during merge; non-empty requirements remain.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,\nS:E484K,escape,S:K417N\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "missing_link",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K",
                "aaDeletions": "",
                "aaInsertions": "",
            },
            {
                "seqName": "has_link",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N",
                "aaDeletions": "",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()
    assert len(reports[1].hits) == 1
    assert reports[1].hits[0].mutation == "S:E484K"


def test_range_autolink_and_user_linked_mutations_compose(tmp_path: Path) -> None:
    """
    Range auto-linking (full range required) AND user linked_mutations should both apply.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:del87-89,range_hit,S:E484K\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            # Full range present but missing linked substitution -> should be gated off
            {
                "seqName": "full_range_missing_sub",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "",
                "aaDeletions": "S:N87-,S:N88-,S:N89-",
                "aaInsertions": "",
            },
            # Full range present and linked substitution present -> should match
            {
                "seqName": "full_range_with_sub",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K",
                "aaDeletions": "S:N87-,S:N88-,S:N89-",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()

    muts = {h.mutation for h in reports[1].hits}
    assert muts == {"S:87-", "S:88-", "S:89-"}


# -----------------------------------------------------------------------------
# Gotcha tests: deletion range edge cases
# -----------------------------------------------------------------------------


def test_deletion_range_start_equals_stop(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del87-87,singleton\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    assert set(d.keys()) == {"S:87-"}
    # When full-range mode is on, linked_mutations should just be the singleton
    assert d["S:87-"].get("linked_mutations", "") == "S:87-"


def test_deletion_range_non_s_gene_expands(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nORF1a:del100-102,orf_hit\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    assert "ORF1a:100-" in d
    assert "ORF1a:101-" in d
    assert "ORF1a:102-" in d


def test_deletion_range_largeish_expansion(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del1-50,big\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    # 50 expanded keys
    expanded = [k for k in d.keys() if k.startswith("S:") and k.endswith("-")]
    assert len(expanded) == 50
    assert "S:1-" in d and "S:50-" in d


def test_deletion_range_permissive_format_allows_spaces(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del 87 - 89,range_hit\n",
    )

    d = load_annotation_table(anno, delimiter=",")

    # Expanded keys exist
    assert "S:87-" in d
    assert "S:88-" in d
    assert "S:89-" in d

    # Auto-linking (default) still requires full range and uses ';' serialization
    assert d["S:87-"].get("linked_mutations") == "S:87-;S:88-;S:89-"
    assert d["S:88-"].get("linked_mutations") == "S:87-;S:88-;S:89-"
    assert d["S:89-"].get("linked_mutations") == "S:87-;S:88-;S:89-"

    # Confirm annotation copied
    assert d["S:87-"]["phenotype"] == "range_hit"
    assert d["S:88-"]["phenotype"] == "range_hit"
    assert d["S:89-"]["phenotype"] == "range_hit"


def test_deletion_range_rejects_stop_less_than_start(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del 89 - 87,bad\n",  # matches permissive regex; invalid ordering
    )
    with pytest.raises(Exception):
        load_annotation_table(anno, delimiter=",")


def test_deletion_range_rejects_zero_start(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del0-2,bad\n",  # matches regex; start <= 0
    )
    with pytest.raises(Exception):
        load_annotation_table(anno, delimiter=",")


def test_deletion_range_rejects_zero_stop(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del1-0,bad\n",  # matches regex; stop <= 0
    )
    with pytest.raises(Exception):
        load_annotation_table(anno, delimiter=",")


def test_deletion_range_rejects_both_zero(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del0-0,bad\n",  # matches regex; both invalid
    )
    with pytest.raises(Exception):
        load_annotation_table(anno, delimiter=",")


def test_deletion_range_rejects_negative_indices_if_present(tmp_path: Path) -> None:
    """
    This is mostly documentary: the current regex does not match negative numbers,
    so this should NOT be treated as a range key and therefore should not raise
    (unless we later add stricter validation).

    If we want negatives to error, you must add an explicit check before regex match.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del-1-3,bad\n",  # does NOT match (?P<start>\\d+)
    )

    # Current behavior: treated as literal mutation key, so loader should succeed.
    d = load_annotation_table(anno, delimiter=",")
    assert "S:del-1-3" in d


def test_deletion_range_rejects_non_numeric_indices_if_strictness_added_later(
    tmp_path: Path,
) -> None:
    """
    Documentary: non-numeric indices do not match the regex, so currently they
    are treated as literal mutation keys and do not raise.

    If we later want to error on these, add a stricter pre-check.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:delA-B,bad\n",  # does NOT match regex
    )

    d = load_annotation_table(anno, delimiter=",")
    assert "S:delA-B" in d


# -----------------------------------------------------------------------------
# Gotcha tests: stop codons + linked gating uses the mapped universe
# -----------------------------------------------------------------------------


def test_stop_codon_match_with_linked_mutations(tmp_path: Path) -> None:
    """
    Annotation may specify stop as gene:{aa}{pos}* while Nextclade stop codons are gene:pos.
    Ensure matching works and linked gating applies.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:N87*,stop_hit,S:E484K\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "stop_without_link",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "S:87",
                "aaSubstitutions": "",
                "aaDeletions": "",
                "aaInsertions": "",
            },
            {
                "seqName": "stop_with_link",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "S:87",
                "aaSubstitutions": "S:E484K",
                "aaDeletions": "",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()
    assert len(reports[1].hits) == 1
    assert reports[1].hits[0].mutation in {"S:N87*", "S:87"}  # depends on which anno key is stored
    assert reports[1].hits[0].annotations.get("phenotype") == "stop_hit"


# -----------------------------------------------------------------------------
# Gotcha tests: delimiter wiring for TSV annotations
# -----------------------------------------------------------------------------


def test_tsv_annotations_with_linked_and_ranges(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.tsv"
    _write_text(
        anno,
        "mutation\tphenotype\tlinked_mutations\nS:del10-12\trange\tS:E484K\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "no_sub",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "",
                "aaDeletions": "S:N10-,S:N11-,S:N12-",
                "aaInsertions": "",
            },
            {
                "seqName": "with_sub",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K",
                "aaDeletions": "S:N10-,S:N11-,S:N12-",
                "aaInsertions": "",
            },
        ],
    )

    reports = compare_nextclade_to_annotations(
        nxt,
        anno,
        delimiter="\t",
    )
    assert reports[0].hits == tuple()
    muts = {h.mutation for h in reports[1].hits}
    assert muts == {"S:10-", "S:11-", "S:12-"}


def _tokenize_serialized_linked(s: str) -> set[str]:
    """Helper to tokenize serialized linked_mutations (always semicolon-joined)."""
    return {t for t in s.split(";") if t}


# LEGACY - NOT DOCUMENTED BEHAVIOUR EXCEPT IN THESE TESTS
def test_split_linked_mutations_accepts_comma_and_semicolon(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        'mutation,phenotype,linked_mutations\nS:E484K,escape,"S:K417N,S:N501Y;S:T478K"\n',
    )

    d = load_annotation_table(anno, delimiter=",")
    lm = d["S:E484K"]["linked_mutations"]

    # Must serialize deterministically with ';'
    assert lm == "S:K417N;S:N501Y;S:T478K"


def test_split_linked_mutations_trims_whitespace_and_ignores_empty(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        'mutation,phenotype,linked_mutations\nS:E484K,escape," ;  S:K417N  ,,  S:N501Y ; "\n',
    )

    d = load_annotation_table(anno, delimiter=",")
    lm = d["S:E484K"]["linked_mutations"]

    assert lm == "S:K417N;S:N501Y"


def test_merge_linked_mutations_is_sorted_and_deduplicated(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\n"
        "S:E484K,escape,S:N501Y\n"
        "S:E484K,escape,S:K417N\n"
        "S:E484K,escape,S:N501Y\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    lm = d["S:E484K"]["linked_mutations"]

    # Must be sorted lexicographically and unique
    assert lm == "S:K417N;S:N501Y"


def test_duplicate_rows_with_empty_linked_mutations_do_not_remove_requirements(
    tmp_path: Path,
) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,\nS:E484K,escape,S:K417N\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    assert d["S:E484K"]["linked_mutations"] == "S:K417N"


def test_linked_mutations_case_sensitive_exact_matching(tmp_path: Path) -> None:
    """
    Parser is case-agnostic; matching is exact-string based.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:E484K,escape,s:k417n\n",
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "has_uppercase",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "S:E484K,S:K417N",
                "aaDeletions": "",
                "aaInsertions": "",
            }
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()


def test_linked_mutations_deletion_format_mismatch(tmp_path: Path) -> None:
    """
    Linked deletions must use canonical {gene}:{pos}- form.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:87-,del,S:N88-\n",  # non-canonical linked token
    )

    nxt = tmp_path / "nextclade.tsv"
    _make_nextclade_tsv(
        nxt,
        rows=[
            {
                "seqName": "has_both_raw",
                "qc.overallStatus": "good",
                "qc.stopCodons.stopCodons": "",
                "aaSubstitutions": "",
                "aaDeletions": "S:N87-,S:N88-",
                "aaInsertions": "",
            }
        ],
    )

    reports = compare_nextclade_to_annotations(nxt, anno, delimiter=",")
    assert reports[0].hits == tuple()


# -----------------------------------------------------------------------------
# Deletion range edge-case correctness
# -----------------------------------------------------------------------------


def test_deletion_range_start_equals_stop_serializes_cleanly(tmp_path: Path) -> None:
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del87-87,single\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    assert set(d.keys()) == {"S:87-"}
    assert d["S:87-"]["linked_mutations"] == "S:87-"


def test_deletion_range_autolink_serialization_uses_semicolon(tmp_path: Path) -> None:
    """
    Auto-linked ranges must serialize using ';' (never comma).
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype\nS:del10-12,range\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    lm = d["S:10-"]["linked_mutations"]

    assert lm == "S:10-;S:11-;S:12-"
    assert "," not in lm


def test_range_and_manual_link_compose_as_union(tmp_path: Path) -> None:
    """
    Range auto-link + user-specified linked_mutations should union.
    """
    anno = tmp_path / "annotations.csv"
    _write_text(
        anno,
        "mutation,phenotype,linked_mutations\nS:del10-12,range,S:E484K\n",
    )

    d = load_annotation_table(anno, delimiter=",")
    lm = d["S:10-"]["linked_mutations"]

    toks = _tokenize_serialized_linked(lm)
    assert toks == {"S:10-", "S:11-", "S:12-", "S:E484K"}
