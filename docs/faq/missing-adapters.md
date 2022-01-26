---
layout: default
parent: FAQ
title: Missing adapters
---

## What are missing adapter tags `ma` and `ac`?

The `ma` and `ac` tags indicate whether the molecule that produces a CCS
read is missing a SMRTbell adapter on its left/start or right/end. The tags
produced by _ccs_ v6.3.0 and newer are based on the `ADAPTER_BEFORE_BAD` and
`ADAPTER_AFTER_BAD` information in the subread `cx` tag.

### Tag `ac`

Array containing four counts, in order:
- *detected* adapters on left/start
- *missing* adapters on left/start
- *detected* adapters on right/end
- *missing* adapter on right/end

### Tag `ma`

Bitmask storing if an adapter is missing on either side of the
molecule. A value of 0 indicates neither end has a confirmed
missing adapter.
- `0x1` if adapter is missing on left/start
- `0x2` if adapter is missing on right/end
