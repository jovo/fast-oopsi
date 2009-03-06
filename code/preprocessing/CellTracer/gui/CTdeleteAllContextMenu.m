function CTdeleteAllContextMenu()
h = gcf;
handles = guidata(h);
hc = findall(h,'Type','uicontextmenu')
if ~isempty(hc)
    for i = 1: length(hc)
        delete(hc(i))
    end
end
hc = findall(h,'Type','uicontextmenu')
%guidata(h,handles);


