

%%
clear;

ranks = [10,10,10,10,20,40,60,80,100,120,140,160,180,200];
tag = [];
algos = {'CPQR','CPQR2pass','LUPP','LUPP2pass','RSVDDEIM','RSVDLS'};
test_CUR_large(1e4, ranks, algos, tag)

%%