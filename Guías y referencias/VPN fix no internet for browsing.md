---
tags:
  - how-to
---

[networking - OpenVPN connecting but no internet access on Ubuntu 16.04 / 18.04 / 20.04 - Ask Ubuntu](https://askubuntu.com/questions/655806/openvpn-connecting-but-no-internet-access-on-ubuntu-16-04-18-04-20-04)

To fix this, edit the OpenVPN connection configuration on Network Manager and click in `IPv4 Settings` tab, then click in `Routes` button.

Then mark `Use this connection only for resources on its network`.

Click `Ok`, then `Save` and reconnect.
