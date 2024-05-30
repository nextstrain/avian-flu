Verify behavior of `transform-field-names`

If the `--field-map` includes a old field name that is the same as the new field
name, this should be no-op.

  $ echo '{"strain": "A"}' \
  >  | $TESTDIR/../../transform-field-names \
  >      --field-map "strain=strain"
  {"strain":"A"}

If the `--field-map` overwrites an existing field, then skip renaming and
print a loud warning.

  $ echo '{"strain": "A", "isolate": "B"}' \
  >  | $TESTDIR/../../transform-field-names \
  >      --field-map "isolate=strain"
  WARNING: skipping rename of isolate because record already has a field named strain.
  {"strain":"A","isolate":"B"}

The `--field-map` may overwrite an existing field if using `--force` flag.

  $ echo '{"strain": "A", "isolate": "B"}' \
  >  | $TESTDIR/../../transform-field-names \
  >      --field-map "isolate=strain" \
  >      --force
  {"strain":"B"}
