<p align="center">
  <img src="doc/img/unanimity.png" alt="unanimity logo"/>
</p>
<h1 align="center">Unanimity</h1>
<p align="center">C++ library and its applications to generate and process accurate consensus sequences</p>

***
## Availability
Latest `ccs` can be installed via bioconda package `pbccs`.

Please refer to our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda)
for information on Installation, Support, License, Copyright, and Disclaimer.

## [Circular Consensus Calling](doc/PBCCS.md)

`ccs` takes multiple reads of the same SMRTbell sequence and combines
them, employing a statistical model, to produce one high quality consensus sequence.
More information available [here](doc/PBCCS.md).

## FAQ

### [Help! I am getting "Unsupported chemistries found: (...)"!](#model-data)

Similar to [pbbam](https://github.com/PacificBiosciences/pbbam), Unanimity's consensus
models require chemistry-dependent parameters. As part of ongoing development efforts,
we might need to introduce new part numbers to identify novel reagents and/or SMRT Cells.
If your version of Unanimity significantly predates the chemistry you have used for
generating collections of data, you will run into issues with Unanimity not being able
to handle your data. In such cases, download the latest version of the model parameters
and place them in a subdirectory of `${SMRT_CHEMISTRY_BUNDLE_DIR}`:

  ```sh
  cd <some persistent dir, preferably the same as used for pbbam>
  export SMRT_CHEMISTRY_BUNDLE_DIR="${PWD}"

  mkdir -p arrow
  cp /some/download/dir/model.json arrow/
  ```

This will cause Unanimity to try to load models from all files in `${SMRT_CHEMISTRY_BUNDLE_DIR}/arrow`
with a `.json` suffix.

## License
[PacBio open source license](LICENSE)

DISCLAIMER
----------
THIS WEBSITE AND CONTENT AND ALL SITE-RELATED SERVICES, INCLUDING ANY DATA, ARE PROVIDED "AS IS," WITH ALL FAULTS, WITH NO REPRESENTATIONS OR WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED, INCLUDING, BUT NOT LIMITED TO, ANY WARRANTIES OF MERCHANTABILITY, SATISFACTORY QUALITY, NON-INFRINGEMENT OR FITNESS FOR A PARTICULAR PURPOSE. YOU ASSUME TOTAL RESPONSIBILITY AND RISK FOR YOUR USE OF THIS SITE, ALL SITE-RELATED SERVICES, AND ANY THIRD PARTY WEBSITES OR APPLICATIONS. NO ORAL OR WRITTEN INFORMATION OR ADVICE SHALL CREATE A WARRANTY OF ANY KIND. ANY REFERENCES TO SPECIFIC PRODUCTS OR SERVICES ON THE WEBSITES DO NOT CONSTITUTE OR IMPLY A RECOMMENDATION OR ENDORSEMENT BY PACIFIC BIOSCIENCES.
