function mod_graph(xp,zp,cizro,alp1,xelek,iter,misfit,nel,ax)
        %         axes(handles.ModRes);
        axes(ax);
        
        deffont='Microsoft Sans Serif';
        hh=patch(xp',zp',repmat(cizro,4,1),'tag','model');
%         handles.model=hh;
        set(hh,'EdgeColor',[170 170 170]/255)
        
        %                 set(hh,'EdgeAlpha',.8)
        alpha(hh,alp1)
        z=unique(zp(:));
        
        for k=1:length(z)
            ZL{k}=sprintf('%6.1f',z(k));
        end
        set(gca,'XLim',[xelek(1)-.01 xelek(end)],'XTick',xelek(1:2:end))
        set(gca,'YLim',[min(zp(:))-.01 max(zp(:))],'YTick',sort(z),'YTickLabel',(ZL))
        if length(z)>8
            set(gca,'FontSize',9)
        else
            set(gca,'FontSize',11)
        end
        
        hpa=colorbar('peer',gca,'Location','eastoutside','FontName',deffont,'tag','colorbar');
%         handles.hpa=hpa;
        xlabel('Distance (m)','FontName',deffont,'FontSize',11);
        ylabel('Depth (m)','FontName',deffont,'FontSize',11);
        title(['Model Resistivity Section  ','Iteration : ',num2str(iter),' RMS = ',sprintf('%5.2f',misfit),' %'] ,'FontName',deffont,'FontSize',11);
%         if ax==handles.ModRes
%             set(handles.ModRes,'Visible','on');
%         end
        
%         imagemenu_tr_patch(handles);
        h(1)=xlabel('Distance (m)','FontName',deffont,'FontSize',11);
        h(2)=ylabel('Depth (m)','FontName',deffont,'FontSize',11);
        h(3)=title(['Model Resistivity Section  ','Iteration : ',num2str(iter),' RMS = ',sprintf('%5.2f',misfit),' %'] ,'FontName',deffont,'FontSize',11);
%         resizeFcn;
        
        
        
    end