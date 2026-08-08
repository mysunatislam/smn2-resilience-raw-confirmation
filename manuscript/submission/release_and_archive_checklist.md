# Permanent Release and Archive Checklist

Do not create the submission tag until the manuscript is genuinely frozen.
The current package is blocked by missing experimental evidence or an explicit
decision to submit without it, and by missing author metadata.

1. Resolve and document the functional-validation and RT-PCR outcome without
   changing the frozen candidate or splice-event framework.
2. Insert author-approved names, affiliations, corresponding-author email,
   ORCIDs, funding, competing interests, contributions, ethics determination,
   and acknowledgements.
3. Regenerate the final manuscript and supplementary PDF, rerun every test, and
   verify the submission manifest hashes.
4. Connect the GitHub repository to Zenodo before publishing the release. The
   official workflow is documented at
   <https://help.zenodo.org/docs/github/enable-repository/>.
5. Create an annotated tag such as `v1.0.0-submission` at the validated commit.
6. Publish a GitHub release from that tag and attach the final supplementary
   PDF plus the machine-readable table archive. GitHub releases are tied to Git
   tags: <https://docs.github.com/en/repositories/releasing-projects-on-github/about-releases>.
7. Confirm that Zenodo archived the release and issued both version-specific
   and concept DOI records.
8. Add the version-specific DOI and tagged GitHub URL to Data and Code
   Availability, the cover letter, and the submission metadata.
9. Rerun the numerical and structure audits after inserting the DOI; do not
   move the tag after publication.
