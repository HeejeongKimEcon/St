function fixeps(h,filename)
% Exports a figure to an eps file while cropping (most of) the white space 
% around. It is also useful when working on several computers, as exported 
% figures are identical across machines.
% Usage:
% fixeps(h,filename)
% Inputs:
% h        - Figure hande
% filename - Name for pdf file
% 
% Christian Bustamante
% October 15, 2021
%
set(h,'renderer','Painters');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,filename,'-depsc','-r0')