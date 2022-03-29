function inp1 = TSS_weekly_input_struct(inp,t)

inp_TS = TSS_input_week(t); % only week t
N = inp.N;

inp1 = inp;
inp1.storeindex = inp_TS.storeindex;
inp1.price = inp_TS.price;
inp1.x = inp_TS.x;
inp1.xp = inp_TS.xp;
inp1.z = inp_TS.z;
inp1.chain0 = squeeze(inp_TS.chain(:,:,1,:));
inp1.T = 1;
inp1.NT = inp.N;
inp1.nu3 = inp.nu3_78;
inp1.on = ones(N,1);
inp1.chain = inp_TS.chain;
inp1.chain2 = inp_TS.chain2;
inp1.ACi2 = inp_TS.ACi2;