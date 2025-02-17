function [] = Plot_GEO(GEO)

%% Plot_GEO function plots the Geometry of a cascaded cylindrical waveguide structure or a cascaded conical waveguide structure
% Inputs: 
%       GEO:       It is a structure input having the following fileds
%       GEO.type - 1 for cascaded cylindrical waveguide structure, 2 for
%                  cascaded conical waveguide structure with a cylindrical waveguide
%                  as the base
%       GEO.N    - If GEO.type is 1, it is the number of the cylindrical
%                  elements. 
%                  If GEO.type is 2, it is the number of conical
%                  elements. So for a conical waveguide structure, the
%                  number of array elements in GEO.E() is N + 1 as it also
%                  has to contain the dimensions of the base waveguide.
%       GEO.E() -  Structure Array with the following fields:
%               -  GEO.E().r - Radius of every element
%               -  GEO.E().l - Length of every element 
% Outputs: 
%        Doesn't return anything. 
%        Only plots the waveguide structure in a surface 2D plot.
% For a proper usage, please check REF.m in the parent directory of the
% tool
%% 
        if GEO.type == 1
            disp('=======================================================');
            disp('Plotting the Antenna Geometry');
            for i = 1:GEO.N
                r(i) = GEO.E(i).r;
                l(i) = GEO.E(i).l;
                [x(i , :, :), y(i, :, :), z(i, :, :)] = cylinder(r(i));
                if i == 1
                    figure; hold on; surf(squeeze(x(i, :, :)), squeeze(y(i, :, :)), squeeze(z(i, :, :))*l(i))
                elseif (i ~= GEO.N) || (i == 2)
                    figure(1); hold on; surf(squeeze(x(i, :, :)), squeeze(y(i, :, :)), sum(l(1:i-1)) + squeeze(z(i, :, :))*l(i))
                    
                    
                    [r1,~] = min([r(i - 1) r(i)]);
                    [r2,~] = max([r(i - 1) r(i)]);
                                   
                    [xc, yc, zc] = cylinder([r1 r2]);
                    
                    hold on; surf(xc, yc, ones(size(zc)) .* sum(l(1:i-1)));
                    
                elseif (i == GEO.N) && (i > 2)
                    hold on; surf(squeeze(x(i, :, :)), squeeze(y(i, :, :)), sum(l(1:i-1)) + squeeze(z(i, :, :))*l(i))
                    
                    [r1,~] = min([r(i - 1) r(i)]);
                    [r2,~] = max([r(i - 1) r(i)]);
                                   
                    [xc, yc, zc] = cylinder([r1 r2]);
                    
                    hold on; surf(xc, yc, ones(size(zc)) .* sum(l(1:i-1)));
                end
            end
           colormap('copper'); % shading('interp')
        else
            for i = 1:GEO.N+1
                r(i) = GEO.E(i).r;
                l(i) = GEO.E(i).l;
                
                if i == 1
                    [x(i , :, :), y(i, :, :), z(i, :, :)] = cylinder(r(i));
                    figure; hold on; surf(squeeze(x(i, :, :)), squeeze(y(i, :, :)), squeeze(z(i, :, :))*l(i))
                else
                    [x(i , :, :), y(i, :, :), z(i, :, :)] = cylinder([r(i - 1) r(i)]);
                    
                    hold on; surf(squeeze(x(i, :, :)), squeeze(y(i, :, :)), sum(l(1:i-1)) + squeeze(z(i, :, :))*l(i))
                    
                end
            end
        end
        colormap('copper'); % shading('interp')
        
        xlabel('X [m]', 'FontSize', 14, 'FontWeight', 'bold');
        ylabel('Y [m]', 'FontSize', 14, 'FontWeight', 'bold');
        zlabel('Z [m]', 'FontSize', 14, 'FontWeight', 'bold');
end
