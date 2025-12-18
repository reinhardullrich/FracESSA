# Prompt to Regenerate Codebase Files

Use this exact prompt to regenerate `codebase_cpp.txt` and `codebase_python.txt`:

```
redo the codebase process, but make 2 files:
1. _cpp: take everything from folder cpp except folder tests, and include cmakelists and release.yml
2. _python: like before, everything from the python folder, but not the results and not the basline .csv file!

3. make a prompt i can run again to get the exact same result!
```

Or simply run the script:

```bash
./analysis/generate_codebase.sh
```

