Setting up SSH key-based authentication allows to avoid entering the password when connecting to it, including in VSCode.

Generate a key pair in the local system:

```
ssh-keygen -t ed25519 -b 4096
```

This would generate two files under the names: `id_ed25519` and `id_ed25519.pub` in the `~/.ssh/` directory.

Copy the key:

```
ssh-copy-id -i ~/.ssh/id_ed25519.pub user@server_ip
```

## Explanation

This command generates a new SSH key pair using the Ed25519 algorithm, which is currently considered secure and efficient. The `-b` flag specifies the number of bits in the key, and 4096 bits is a common choice for key length.

**Key Pair Usage:**

- The private key (`id_ed25519`) stays on your local machine.
- The public key (`id_ed25519.pub`) is copied to the server.

When you attempt to connect to the server, the server checks if the private key on your local machine matches any of the public keys stored in the `authorized_keys` file on the server. If there's a match, you're granted access without entering a password.
## References

- [Visual Studio Code Remote Development Troubleshooting Tips and Tricks](https://code.visualstudio.com/docs/remote/troubleshooting#_ssh-tips)
- [visual studio code - How to save ssh password to vscode? - Stack Overflow](https://stackoverflow.com/questions/66113731/how-to-save-ssh-password-to-vscode)
