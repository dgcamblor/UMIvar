
To make Markdown files in a Git repository detectable by Linguist for language statistics, you can use a `.gitattributes` file to override the default behavior. Here's how you can do it:

1. Create a file named `.gitattributes` in the root of your repository if it doesn't exist already.
2. Add the following line to the `.gitattributes` file: `*.md linguist-detectable=true`.

There are more settings that can be found at: [linguist/docs/overrides.md at master · github-linguist/linguist · GitHub](https://github.com/github-linguist/linguist/blob/master/docs/overrides.md#using-gitattributes).