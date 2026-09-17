# Maintained Analysis Code

This directory contains the repository's root-level analysis scripts. Run scripts from the repository root when they use project-relative inputs, for example:

```bash
python src/count_ctcf_in_domains.py
```

New code should take input and output paths as explicit arguments or configuration values rather than relying on the current working directory.