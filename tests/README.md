# Tests

The following tests with output are included:

- 1. Labelfree
- 1a. Labelfree with warnings
- 2. Labelfree / phospho
- 2a Labelfree / phospho / add a set
- 2b. Labelfree / phospho / no decoy
- 2c. Labelfree / phospho / no target
- 3. TMT16
- 3a TMT16 / add a set
- 3b TMT16 mixed fractionated / non fractionated
- 3c TMT16 single file (tests quoting alright)
- 3d TMT16 / phospho+acetyl / PTM total protein normalization
- 4. TMT18 / phospho / median sweep
- 4a TMT18 / phospho / median sweep / add a set
- 4b TMT18 / phospho / median sweep / rerun with different post-PSM parameters
- 3+4a TMT16/18 mixed isobaric types
- 5. TIMS instrument (different MS1 centroid)

These tests will output with an expected error:
- 2c. Labelfree / phospho / no target

Some tests depend on eachother in order to not have to rerun certain expensive steps. In several tests the filenames, setnames etc will contain characters that should be quoted or escaped.

## TODO:
- github actions
- error checking
- output checking
- test that crashes TMT / DEqMS errors
