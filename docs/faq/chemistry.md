---
layout: default
parent: FAQ
title: Chemistry
---

## Help! I am getting "Unsupported ..."!
If you encounter the error `Unsupported chemistries found: (...)` or
`unsupported sequencing chemistry combination`, your _ccs_ binaries do not
support the used sequencing chemistry kit, from here on referred to as
"chemistry". This may be because we removed support of an older chemistry or
your binary predates release of the used chemistry. This is unlikely to happen
with _ccs_ from SMRT Link installations, as SMRT Link is able to automatically
update and install new chemistries. Thus, the easiest solution is to always use
_ccs_ from the SMRT Link version that shipped with the release of the sequencing
chemistry kit.

**Old chemistries:** With _ccs_ 4.0.0, we have removed support for the last RSII
chemistry `P6-C4`. The only option is to downgrade _ccs_ with `conda install
pbccs==3.4`.

**New chemistries:** It might happen that your _ccs_ version predates the
sequencing chemistry kit. To fix this, install the latest version of _ccs_ with
`conda update --all`. If you are an early access user, follow the [monkey patch
tutorial](/faq/chemistry#monkey-patch-ccs-to-support-additional-sequencing-chemistry-kits).

## Monkey patch _ccs_ to support additional sequencing chemistry kits
Please create a directory that is used to inject new chemistry information
into _ccs_:

```sh
mkdir -p /path/to/persistent/dir/
cd /path/to/persistent/dir/
export SMRT_CHEMISTRY_BUNDLE_DIR="${PWD}"
mkdir -p arrow
```

Execute the following step by step instructions to fix the error you are observing
and afterwards proceed using _ccs_ as you would normally do. Additional chemistry
information is automatically loaded from the `${SMRT_CHEMISTRY_BUNDLE_DIR}`
environmental variable.

### Error: "unsupported sequencing chemistry combination"
Please download the latest out-of-band `chemistry.xml`:

```sh
wget https://raw.githubusercontent.com/PacificBiosciences/pbcore/develop/pbcore/chemistry/resources/mapping.xml -O "${SMRT_CHEMISTRY_BUNDLE_DIR}"/chemistry.xml
```

### Error: "Unsupported chemistries found: (...)"
Please get the latest consensus model `.json` from PacBio and
copy it to:

```sh
cp /some/download/dir/model.json "${SMRT_CHEMISTRY_BUNDLE_DIR}"/arrow/
```
