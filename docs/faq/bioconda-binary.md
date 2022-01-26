---
layout: default
parent: FAQ
title: Bioconda binary
---

## The binary does not work on my linux system!
Contrary to official SMRT Link releases, the `ccs` binary distributed via bioconda
is tuned for performance while sacrificing backward compatibility.
We are aware of following errors and limitations. If yours is not listed, please
file an issue on our [official pbbioconda page](https://github.com/PacificBiosciences/pbbioconda).

**`Illegal instruction`** Your CPU is not supported.
A modern (post-2008) CPU with support for
[SSE4.1 instructions](https://en.wikipedia.org/wiki/SSE4#SSE4.1) is required.
SMRT Link also has this requirement.

**`FATAL: kernel too old`** Your OS or rather your kernel version is not supported.
Since _ccs_ v4.2 we also ship a second binary via bioconda `ccs-alt`, which does
not bundle a newer `glibc`. Please use this alternative binary.

For _ccs_, we offer two binaries in bioconda:

 * `ccs`, statically links `glibc` v2.33 and `mimalloc` v1.6.3.
 * `ccs-alt`, was build by dynamically linking `glibc` v2.17, but statically links `mimalloc` v1.6.3.
