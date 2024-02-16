function fixpdf(h,filename)
% Exports a figure to a pdf file while cropping (most of) the white space 
% around. It is also useful when working on several computers, as exported 
% figures are identical across machines.
% Usage:
% fixpdf(h,filename)
% Inputs:
% h        - Figure hande
% filename - Name for pdf file
% 
% Christian Bustamante
% November 12, 2018
%
set(h,'renderer','Painters');
set(h,'Units','Inches');
pos = get(h,'Position');
set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])
print(h,filename,'-dpdf','-r0')