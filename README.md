# dynamic_control_wlan
Implementation of "Dynamic Control for On-demand Interference-managed WLAN Infrastructures"

# Prerequisites
- JAVA
- CPLEX solver
- JSON (JavaScript Object Notation)
- Odin Framework
-- Odin client which inclues click should be installed in each AP

# Codes
- JITWLAN.java: implementation of main algorithm in the paper
- Test_JIT_WLAN.java: sample test code
- odin-master.7z : Whole JITWLAN code which is included in odin-master

# Start up setting
-Install Odin-framework
See https://github.com/lalithsuresh/odin
- Add-client : Add client's mac-address and IP to client file
(In this file abc.txt, original Odin using "net.floodlightcontroller.odin.master.OdinMaster.clientList")
- Add-switch : Add switch's address to switch file
(In this file zxc.txt, original Odin using "net.floodlightcontroller.odin.master.OdinMaster.poolFile")

-Basic AP setting



