---
tags:
  - how-to
---

- Download the [ME3 coalesced utility](https://wenchy.net/me3-coalesced-utility/) 
- Follow the instructions in the reddit post to find the Tab binding and remap PC_Push_to_Talk to something else (e.g., F10)
- Use the utility's Find feature to find all references to "binding" and then every time "B" or "T" is listed, remap them to something else (e.g., F11 and F12).

For anyone else with this issue: I had to download an editor (https://wenchy.net/me3-coalesced-utility/) for the Coalesced.bin file, which is found in the ME3/BioGame/CookedPCConsole directory.

Then I had to navigate to the right settings, bioinput.ini > sfxgame > sfxgamemodedefault > bindings.

Then in the entry there, I had to change the PC_Push_to_Talk entry from Tab to F10, and the PC_Command_Console to Tab.