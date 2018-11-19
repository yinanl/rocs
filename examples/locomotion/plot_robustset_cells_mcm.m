function plot_robustset_cells_mcm(N1,N2,xpos,ds,dz,omega,x,z0,zmin,s)
% FUNCTION
% -----------
% plot the partition of a robust margin set for MCM mode.
%
% INPUT
% -----------
% N1,N2: the number of grids along sigma and zeta, respectively.
% xpos: an (N x 1) array of x values in Euclidean space.
% ds, dz: the granularity for the partition along sigma and zeta,
%         respectively.
% omega: a fixed control value.
% x: the keyframe state.
% z0: zeta value for the keyframe state.
% zmin: the lower bound of zeta for all cells.
% s: deceleration ('left') or acceleration ('right') of MCM.


global CP LThin LThick;

switch s
    case 'left'
        coeff= log(2);
        dz0 = dz;
    case 'right'
        coeff= 1/log(2);
        dz0= log(2)^(N2-1)*dz;
end

for i= 1:N1
    dz= dz0;
    for j= 1:N2
        dsl= (-N1+2*(i-1))*ds;
        dsu= (-N1+2*(i-1) + 2)*ds;
        if (j<2)
            dzl= zmin;
            dzu= dzl + dz;
            dz= coeff*dz;
        else
            dzl= dzu;
            dzu= dzl + dz;
            dz= coeff*dz;
        end
        
        xtipl= dsl/(2*omega) + x(1);
        xtipu= dsu/(2*omega) + x(1);
        % sigma tube
        vl= sqrt(abs(2*omega*(xpos-x(1))-dsl+x(2)^2));
        vu= sqrt(abs(2*omega*(xpos-x(1))-dsu+x(2)^2));
        % zeta tube
        zl= x(2)*exp((dzl/z0-1 + xpos-x(1))*(-1/omega));
        zu= x(2)*exp((dzu/z0-1 + xpos-x(1))*(-1/omega));
        
        % % plot robust tubes
%         switch s
%             case 'left'
%                 plot(xpos(xpos<=xtipu), vu(xpos<=xtipu), '-',...
%                     'Color', CP(5,:), 'LineWidth', LThin)
%                 plot(xpos(xpos<=xtipl), vl(xpos<=xtipl), '-',...
%                     'Color', CP(5,:), 'LineWidth', LThin)
%             case 'right'
%                 plot(xpos(xpos>=xtipu), vu(xpos>=xtipu), '-',...
%                     'Color', CP(5,:), 'LineWidth', LThin)
%                 plot(xpos(xpos>=xtipl), vl(xpos>=xtipl), '-',...
%                     'Color', CP(5,:), 'LineWidth', LThin)
%             otherwise
%                 
%         end
%         plot(xpos, zl, '-', 'Color', CP(5,:), 'LineWidth', LThin)
%         plot(xpos, zu, '-', 'Color', CP(5,:), 'LineWidth', LThin)
        
        % % plot cells
        zero_vector= zeros(size(xpos,1),1);
        switch s
            case 'left'
                vu_1= [xpos [vu(xpos<=xtipu); zero_vector(xpos>xtipu)]];
                vl_1= [xpos [vl(xpos<=xtipl); zero_vector(xpos>xtipl)]];
                lline= vu_1(vu_1(:,2)>=zl & vu_1(:,2)<=zu,:);
                rline= vl_1(vl_1(:,2)>=zl & vl_1(:,2)<=zu,:);
            case 'right'
                vu_1= [xpos [zero_vector(xpos<xtipu); vu(xpos>=xtipu)]];
                vl_1= [xpos [zero_vector(xpos<xtipl); vl(xpos>=xtipl)]];
                lline= vu_1(vu_1(:,2)>=zu & vu_1(:,2)<=zl, :);
                rline= vl_1(vl_1(:,2)>=zu & vl_1(:,2)<=zl, :);
            otherwise
        end
        tline= [xpos(zu>=vu_1(:,2) & zu<=vl_1(:,2)),...
            zu(zu>=vu_1(:,2) & zu<=vl_1(:,2))];
        bline= [xpos(zu>=vu_1(:,2) & zu<=vl_1(:,2)),...
            zl(zu>=vu_1(:,2) & zu<=vl_1(:,2))];
        plot(lline(:,1),lline(:,2), '-',...
            'Color', CP(2,:),'LineWidth', LThick)
        plot(rline(:,1), rline(:,2), '-',...
            'Color', CP(2,:),'LineWidth', LThick)
        plot(tline(:,1), tline(:,2), '-',...
            'Color', CP(2,:),'LineWidth', LThick)
        plot(bline(:,1), bline(:,2), '-',...
            'Color', CP(2,:),'LineWidth', LThick)
    end
end