Verify behavior of `transform-genbank-location` around prescence/abscence of
`database` and `location` fields.

If `location` field is present, transform it.

  $ echo '{"database":"GenBank", "location": "USA:Oregon, Salem" }' \
  >  | $TESTDIR/../../transform-genbank-location
  {"database":"GenBank","location":"Salem","country":"USA","division":"Oregon"}

If `database` field is absent, complain.

  $ echo '{"location": "USA:Oregon, Salem" }' \
  >  | $TESTDIR/../../transform-genbank-location
  Record must contain `database` field to use `transform-genbank-location.`
  {"location":"USA:Oregon, Salem"}

If `database` field has unsupported value, complain.

  $ echo '{"database": "unsupported", "location": "USA:Oregon, Salem" }' \
  >  | $TESTDIR/../../transform-genbank-location
  Database value of unsupported not supported for `transform-genbank-location`; must be "GenBank" or "RefSeq".
  {"database":"unsupported","location":"USA:Oregon, Salem"}


If `location` field is absent, complain.

  $ echo '{"database": "GenBank" }' \
  >  | $TESTDIR/../../transform-genbank-location
  `transform-genbank-location` requires a `location` field; this record does not have one.
  {"database":"GenBank"}
