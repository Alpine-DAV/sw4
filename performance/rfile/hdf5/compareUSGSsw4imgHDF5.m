% Test for Rfile hdf5 output in SW4
% Assumes you have run the tests in the sw4/tool/essi directory:
%   sw4 berkeley-small-att-h100.in (produces berkeley-small-att-h100/ files)
clear;

% basedir = pwd+"/cori";
basedir = pwd;

% ---------------------------------------------------------------------
% entire-USGS-h400
% ---------------------------------------------------------------------
testname = 'entire-berkeley';
% fig_axis = [10000 12000 -2000 500]; % berkeley zoom
fig_axis = [0 12000 -6387.5 500]; % berkeley entire
% testname = 'entire-usgs';
% fig_axis = [210000 240000 -8000 1000];
fprintf('%%% Beginning %s test\n', testname);

cycle = 0; % cycle at which we'll compare output
cycle_str = sprintf('%01d',cycle); % string cycle number
% yslice = 64000;
% yslice = 6000;
yslice = 0;
yslice_str = sprintf('%d',yslice); % yslice as a string
file = sprintf('%s/%s.sfile', basedir, testname);
if (exist(file))
    h_coarse = h5read(file,'/Coarsest horizontal grid spacing');
    h5_ngrids = h5read(file,'/ngrids');
    h5_gridnz = h5read(file,'/grid nz');
    mat_h5 = cell(h5_ngrids+1,10); % extra for topo, 5 comps + 5 other vars
    
    for g=0:h5_ngrids-1
        loc = sprintf('/Material_model/grid_%d',g);
        grid_id = h5readatt(file,loc,'grid_id');
        assert(grid_id == g);

        comp_str = ["Rho", "Cp", "Cs", "Qp", "Qs"];
        for comp=1:5
            data = h5read(file,sprintf('%s/%s',loc,comp_str(comp)));
            % NB: matlab indices are reversed k,j,i
            mat_h5{g+1,comp} = permute(data,[3 2 1]);
        end
    end
    
    % Read the topo for each grid interface
    for g=0:h5_ngrids
        loc = sprintf('/Z_interfaces/z_values_%d',g);
        interface_id = h5readatt(file,loc,'interface_id');
        assert(interface_id == g);

        data = h5read(file,loc);
        % NB: matlab indices are reversed k,j,i
        mat_h5{g+1,6} = permute(data,[2 1]);
    end
    
    % Read in the appropriate 3Dimg file
    img_ngrids = 1; % Will read it from the file  
    comp_str = ["rho", "p", "s", "qp", "qs"];
    g=1;
    while (g<=img_ngrids)
        for comp=1:5
            file = sprintf('%s/%s-results/image.cycle=%s.y=%s.%s.sw4img', ...
                basedir,testname,cycle_str,yslice_str,comp_str(comp));
            if (~exist(file))
                warning('Unable to open .sw4img file: %s', file); 
                clear mat_img;
                break;
            end
%             [im,x,y,z,t,timestring,zmin]=readimage3d(file,g,0);
            [im,x,y,z,plane,t,timestring,npatches,zmin]=readimage(file, g, 1);
            img_ngrids = npatches;            
            if ~exist('mat_img')
                mat_img = cell(img_ngrids+1,10); % same read from .3Dimg files
            end
            mat_img{g,comp} = im;
            % Save topo too
            if (comp == 1)
                if (size(z,1) == 1) % Cart, no curv
                    mat_img{g+1,6} = (zmin+z(1)) * ones(size(x,2),size(y,2));
                else % Curv
                    mat_img{g+1,6} = z(:,:,1);
                end
                % Save a version of x, y, z full coordinates for plotting
                mat_img{g,7} = x;
                mat_img{g,8} = y;
                mat_img{g,9} = z; % z+zmin;
                
                if (g == 1)
                    % Save bottom interface, too
                    mat_img{g,6} = (zmin + z(size(z,2))) * ones(size(x,2),size(y,2));
                end
            end
        end
        g = g+1;
    end
    
    % plot a vertical slice of each variable
    comp_str = ["Rho", "Cp", "Cs", "Qp", "Qs"];
    for fig=1:5
        figure(fig)
        clf;
    end
    
    % calculate min/max for all comps across grids
    min_vals = 1e38*ones(5,1);
    max_vals = -1e38*ones(5,1);
    for g=1:img_ngrids
        for comp=1:5
            data = mat_img{g,comp};
            min_vals(comp) = min(min_vals(comp), min(data,[],'all'));
            max_vals(comp) = max(max_vals(comp), max(data,[],'all'));
        end
    end

    for g=1:h5_ngrids
        for comp=1:5
            data = mat_h5{g,comp};
            min_vals(comp) = min(min_vals(comp), min(data,[],'all'));
            max_vals(comp) = max(max_vals(comp), max(data,[],'all'));
        end
    end
    
    % plot the img data
    ncontours = 20;
    t = (0:ncontours)/ncontours;
    contour_vals = round(min_vals*(1-t) + max_vals*t);
    x = mat_img{1,7};
    h = x(2)-x(1);
    crs_slicej = round(yslice/h)+1;
    slicej = crs_slicej;
    for g=1:img_ngrids
        zbot = mat_img{g,6};
        ztop = mat_img{g+1,6};
        x = mat_img{g,7};
        h = x(2)-x(1);
        if ((g>1) && (g<img_ngrids))
            h = h/2;
            slicej = (slicej-1)*2+1;
        end
        z = mat_img{g,9};
        nx = size(x,2);
        ny = size(y,2);
        if (size(z,1) == 1)
            % Cartesian
            nz = size(z,2);
            z = ones(nx,1) * z;
        else
            % Topo
            nz = size(z,3);
%             z = reshape(z(:,slicej,:),[nx nz]);
        end
        x = x' * ones(1, nz);
        for comp=1:5
            data = mat_img{g,comp};
%             data = reshape(data(:,slicej,:),[nx nz]);
            figure(comp);
            subplot(1,2,2);
            title(comp_str(comp) + ' sfile');
            subplot(1,2,1);
            title(comp_str(comp) + ' sw4img');
            if (g==img_ngrids)
                zplot = -z';
            else
                zplot = -z;
            end
            contourf(x, zplot, data', contour_vals(comp,:), 'Showtext', 'on');
            axis(fig_axis);
            colorbar;
            grid on;
            hold on;
%             figure(comp+5);
%             subplot(1,2,2);
%             title(comp_str(comp) + ' hdf5');
%             grid on;
%             subplot(1,2,1);
%             title(comp_str(comp) + ' 3Dimg');
%             grid on;
        end
    end
    
    % plot the hdf5 data
    h = h_coarse;
%     crs_slicej = round(yslice/h)+1;
    crs_slicej = yslice/h + 1;
    slicej = round(crs_slicej);
    for g=1:h5_ngrids
        zbot = mat_h5{g,6};
        ztop = mat_h5{g+1,6};
        nx = size(ztop,1);
        ny = size(ztop,2);
        if (g>1)
            h = h / 2;
            crs_slicej = (crs_slicej-1)*2+1;
            slicej = round(crs_slicej);
            % Need to interpolate zbot from coarse interface to this one
            ztmp = zeros(size(ztop));
            cnx = nx/2+1;
            cny = ny/2+1;
            ztmp(1:2:nx,1:2:ny) = zbot; % copy odd points
            % avg even x points
            ztmp(2:2:nx-1,1:2:ny) = .5*(ztmp(1:2:nx-2,1:2:ny)+ztmp(3:2:nx,1:2:ny)); 
            % avg even y points
            ztmp(1:2:nx,2:2:ny-1) = .5*(ztmp(1:2:nx,1:2:ny-2)+ztmp(1:2:nx,3:2:ny)); 
            % avg center points
            ztmp(2:2:nx-1,2:2:ny-1) = .25*(ztmp(1:2:nx-2,1:2:ny-2)+ztmp(1:2:nx-2,3:2:ny) ...
                +ztmp(3:2:nx,1:2:ny-2)+ztmp(3:2:nx,3:2:ny));
            zbot = ztmp;
        end
        nz = cast(h5_gridnz(g),'double');
        x = (0:nx-1)'*h * ones(1,nz);
        % Create z as constant grid spacing in each column between zbot, ztop
        t = (0:nz-1)/(nz-1);
        z = ztop(:,slicej)*(1-t) + zbot(:,slicej)*t;
        for comp=1:5
            data = mat_h5{g,comp};
            data = reshape(data(:,slicej,:),[nx nz]);
            figure(comp);
            subplot(1,2,2);
            title(comp_str(comp) + ' sfile');
            contourf(x, -z, data, contour_vals(comp,:), 'Showtext', 'on');
            colorbar;
            axis(fig_axis);
            grid on;
            hold on;
%             figure(comp+5);
%             subplot(1,2,2);
%             title(comp_str(comp) + ' hdf5');
%             grid on;
%             subplot(1,2,1);
%             title(comp_str(comp) + ' 3Dimg');
%             grid on;
        end
    end
    
    for (f=1:5)
        fh = figure(f);
        fh.Position = [93 67 1538 783];
%         fh = figure(f+5);
%         fh.OuterPosition = [-2180 152 1288 782];
    end    
else
    warning('Skipping %s test, cannot open hdf5 file: %s', testname, file);
end


