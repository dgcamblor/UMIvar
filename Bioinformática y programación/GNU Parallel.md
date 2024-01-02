The most important items are:

- `{}`: Placeholder. It acts as placeholder for the actual values from the input sources.
- `:::`: Input separator. It separates the command from the input sources.

For example:

```
parallel 'command {}' ::: file1 file2 file3
```

Runs a command in several files in parallel.