# hrstrans
General Second Order Transport Matrix Code
Matrix elements are taken from Karl Brown SLAC technote on Beam Optics

To run hrstrans 
Instantiate the class for example, 
THRSTrans *trans = new THRSTrans(0.122396*0.94, -0.136543, -0.171633*0.93, 0.5, 0.5, 0.050178, 0.037056, THRSTrans::kPREX, THRSTrans::kLHRS,evtnumber);

Then you can call trans->tree() to create a root file to do subsequent analysis. 

