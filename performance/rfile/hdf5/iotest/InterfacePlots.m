% Script for plotting sfile interfaces
basedir='.';
testname='berkeley';
file = sprintf('%s/%s.sfile', basedir, testname);
h_coarse = h5read(file,'/Coarsest horizontal grid spacing');
h5_ngrids = h5read(file,'/ngrids');
h5_gridnz = h5read(file,'/grid nz');
mat_h5 = cell(h5_ngrids,2); % 2 comps for top, bot

for b=0:1 % bot/top
for g=0:h5_ngrids-1
    if (b==0)
        field = 'bot';
    else
        field = 'top';
    end
    loc = sprintf('/Z_interfaces/z_%s_%d',field,g);
    attr = sprintf('z_%s_id',field);
    interface_id = h5readatt(file,loc,attr);
    assert(interface_id == g);

    data = h5read(file,loc);
    % NB: matlab indices are reversed k,j,i
    mat_h5{g+1,1+b} = permute(data,[2 1]);
end
end

