# Structural prediction manifest validation

The manifest validator checks structural prediction manifests before submission or retry.

It validates:

- header compatibility
- expected column count
- duplicate job_id values
- empty path columns
- known status values

This is useful because TSV manifests with long output paths are hard to read directly in a terminal. Validation provides a safer check than visual inspection.
