"""docs/SCHEMA.md is generated, so CI has to notice when someone edits the schema without it."""

import generate_schema_doc


def test_schema_doc_is_up_to_date():
    assert generate_schema_doc.main(["--check"]) == 0, (
        "docs/SCHEMA.md is stale; regenerate with "
        "`python rnaseq_helper_scripts/generate_schema_doc.py`"
    )


def test_schema_doc_documents_every_required_table():
    from schemas import REQUIRED_TABLES, table as get_table

    content = generate_schema_doc.DOC_PATH.read_text()
    for name in REQUIRED_TABLES:
        spec = get_table(name)
        assert f"### `{spec.relpath}`" in content
        for column in spec.columns:
            assert f"| `{column.name}` |" in content
