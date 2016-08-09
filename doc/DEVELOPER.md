# Developer environment

## Clang-format pre commit hook

### What?
The intention here is to get source code formatting automatically
checked and guaranteed before checkin.

### Why?
- Automate style-guide compliance to avoid squabbles
- Post-checkin reformatting binges mess up history.  It's best to get
  it right from the get-go.


### How?
See the tools:
  - `tools/check-formatting --staged` for fast style-checking of
    changed files, in a git commit/push hook;
  - `tools/check-formatting --all` for slower style-checking of everything
  - `tools/format-all`to format everything

To enable the checking as a pre-commit hook, place the following in
`.git/hooks/pre-commit` (and make that file executable):

```sh
#!/bin/bash
./tools/check-formatting --staged
```

### Tips
- We should *not* reformat bundled third-party code as it will make it
  difficult to diff with upstream.  Put such code in a subdirectory
  with its own .clang-format file containing: "BasedOnStyle: None" to
  disable formatting.

- If clang-format mangles something (for example, a comment block with
  ASCII art or significant whitespace), you can protect a block using
  `// clang-format off` and then `// clang-format on`.

* What are these binaries?
  The `clang-format` binaries that are checked in are static builds of
  clang-format v3.9 from https://github.com/angular/clang-format.
  They are bundled because the LLVM toolchain has a reputation for
  rapid change and incompatibility.  We want to guarantee that devs
  are using the same version as CI is.
