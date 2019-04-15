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

-Basic AP/User setting : Add user list and switch list at Odinmaster.java

Around line 1181, set "numOfUser", and add UserMacList and switchIPList.

EX) 
userMacList.add(0, MACAddress.valueOf("88:32:9B:9F:54:F1"));
switchIPList.add(0, "192.168.1.51");

Around line 1303
set start up setting number

EX)
int Num_user =8;//cList.size();
int Num_AP =10;
int Num_Channel = 1;
sOnOff = new int[Num_AP];

Around line 1378
setting overlap graph using "solver.setOverlap_G" function
This function has 3 input parameter
first and second is switch's number and third is interference flag.

EX)
solver.setOverlap_G(0, 1, 1);//ap0>ap1 interference
solver.setOverlap_G(1, 0, 1);//ap1>ap0 interference

Around line 1378
set user IP which written at "clientList" and controller IP 

EX)
String serverIP = "143.248.56.10";
String user1IP = "143.248.57.70";


# Start evaluation
1.switch and user add step
After run floodlight, add switch to controller.
Each switch are add to controller in order(this is need for fit overlap graph which write with using solver.setOverlap_G function ).
After switch are added to switch, each user connect to wi-fi which odin_vlan's SSID.
this step must have done in 5min (you can modify this time change number at line 1221 in odinmaster.java file.

2. send data to user
when "Start test. after 0min, flow start!" text is printed in controller consol window, send data to user

3. each 5min period, test is repeat
Our test sinario repeated in 5 min period.
when your test is over(all data sended to user), turn off the controller.
