function plot_robustset_cells_pipm(N1,N2,xpos,ds,dz,omega,x,Dx,x0,z0)
% FUNCTION
% -----------
% plot the partition of a robust margin set for PIPM mode
%
% INPUT
% -----------
% N1,N2: the number of grids along sigma and zeta, respectively.
% xpos: an (N x 1) array of x values in Euclidean space.
% ds, dz: the granularity for the partition along sigma and zeta,
%         respectively.
% omega: a fixed control value.
% x: the keyframe state.
% Dx: Dx = x0-x.
% x0: an initial state different from x.
% z0: zeta value for the keyframe state.


global CP LThin LThick;

for i1= 1:N1
    for j1= 1:N2
        dsl= (-N1+2*(i1-1))*ds;
        dsu= (-N1+2*i1)*ds;
        dzl= (-N2+2*(j1-1))*dz;
        dzu= (-N2+2*j1)*dz;
        
        % sigma tube
        vl= real(sqrt(omega^2*(xpos-x(1)).^2 + x(2)^2 + omega^2/x(2)^2*dsl));
        vu= real(sqrt(omega^2*(xpos-x(1)).^2 + x(2)^2 + omega^2/x(2)^2*dsu));
        % zeta tube
        zl= (dzl/z0-1)*Dx./(xpos-x(1));
        zlp= zl(zl>=0);
        ilp= find(zl>=0);
        xpl= xpos(ilp);
        zlp= x0(2)*zlp.^(1/omega^2);
        zu= (dzu/z0-1)*Dx./(xpos-x(1));
        zup= zu(zu>=0);
        iup= find(zu>=0);
        xpu= xpos(iup);
        zup= x0(2)*zup.^(1/omega^2);
        
        % % plot robust tubes
%         plot(xpos, vu, '-', 'Color', CP(5,:), 'LineWidth', LThin)
%         plot(xpos, vl, '-', 'Color', CP(5,:), 'LineWidth', LThin)
%         plot(xpu, zup, '-', 'Color', CP(5,:), 'LineWidth', LThin)
%         plot(xpl, zlp, '-', 'Color', CP(5,:), 'LineWidth', LThin)
        
        vuu= [xpu vu(iup)];
        vul= [xpl vu(ilp)];
        vlu= [xpu vl(iup)];
        vll= [xpl vl(ilp)];
        zbdu= [xpu zup];
        zbdl= [xpl zlp];
        vbdu= [xpos vu];
        vbdl= [xpos vl];
        % % left and right boundary of the initial robust set
        rbleft= zbdl(vll(:,2)<=zlp & zlp<=vul(:,2),:);
        rbright= zbdu(vlu(:,2)<=zup & zup<=vuu(:,2),:);
        % % upper and lower boundary of the initial robust set
        if (rbleft(end,2)>rbleft(1,2))
            pua= rbleft(end,1);
            pla= rbleft(1,1);
        else
            pua= rbleft(1,1);
            pla= rbleft(end,1);
        end
        if (rbright(end,2)>rbright(1,2))
            pub= rbright(end,1);
            plb= rbright(1,1);
        else
            pub= rbright(1,1);
            plb= rbright(end,1);
        end
        rbupper= vbdu(pua<=vbdu(:,1) & vbdu(:,1)<=pub, :);
        rblower= vbdl(pla<=vbdl(:,1) & vbdl(:,1)<=plb, :);
        
        % % plot the boundary of robust sets
        plot(rbleft(:,1),rbleft(:,2), '-',...
            'Color', CP(2,:),'LineWidth', LThick)
        plot(rbright(:,1), rbright(:,2), '-',...
            'Color', CP(2,:),'LineWidth', LThick)
        plot(rbupper(:,1), rbupper(:,2), '-',...
            'Color', CP(2,:),'LineWidth', LThick)
        plot(rblower(:,1), rblower(:,2), '-',...
            'Color', CP(2,:),'LineWidth', LThick)
    end
end