hAxes=findall(gcf,'type','axes');

if exist('hXlabel','var')
    set(hXlabel , ...
        'FontName'   , 'DejaVuSans');

    set(hXlabel, ...
        'interpreter'   , 'latex', ...
        'FontSize'      , 11          );
end

if exist('hYlabel','var')
    set(hYlabel , ...
        'FontName'   , 'DejaVuSans');

    set(hYlabel , ...
        'interpreter'   , 'latex', ...
        'FontSize'      , 11          );
end



if exist('hColorbar','var')
        set(hColorbar ,...
        'FontName'   , 'DejaVuSans',...
        'FontSize'   , 8);
        hColorbarLabel=get(hColorbar,'YLabel');
        set(hColorbarLabel,'FontName','DejaVuSans');
        set(hColorbarLabel,'FontSize',8);
end


if exist('hLegend','var')
    set(hLegend,'FontName','DejaVuSans','FontSize'   , 11);
    set(hLegend , ...
    'interpreter'   , 'latex');
    hAxes=hAxes(hAxes~=hLegend);
    %set(hLegend,'Location','northeastoutside');
end

if exist('hTitle','var')
set( hTitle                    , ...
    'FontSize'   , 16          , ...
    'FontWeight' , 'normal'      );
set(hTitle , ...
    'FontName'   , 'DejaVuSans');
set(hTitle , ...
    'interpreter'   , 'latex');
end

set( hAxes                       , ...
    'FontName'   , 'DejaVuSans' );


set(hAxes, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'on'      , ...
  'YMinorTick'  , 'on'      , ...
  'LineWidth'   , 1         );
% 
% Adjust linewidth
% hLines = get(hAxes,'children');
% for n=1:length(hLines);
%     if strcmpi(get(hLines(n),'type'),'line')
%         set(hLines(n),'LineWidth',0.75);
%     end
% end
