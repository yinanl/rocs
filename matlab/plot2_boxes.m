function plot2_boxes(box, varargin)
% PLOT2_BOXES(BOX,VARARGIN)
% -------------------------
% Inputs:
% 
if(isempty(varargin) || nargin~=4)
    facecolor= [0 0.4470 0.7410];
    edgecolor= 'k';
    a= 1;
else
    facecolor= varargin{1};
    edgecolor= varargin{2};
    a= varargin{3};
end

if ~isempty(box)
    
    if(isa(box, 'interval'))
        type= 'interval';
    elseif(isa(box, 'numeric'))
        type= 'numeric';
    else
        error('Wrong input type.')
    end
    
    draw(box, type, facecolor, edgecolor, a);
end

alpha(a);

end

function draw(box, type, varargin)

switch type
    case 'interval'
        x= box(:,1);
        y= box(:,2);
        xl= x.lower;
        xu= x.upper;
        yl= y.lower;
        yu= y.upper;
    case 'numeric'
        xl= box(:, 1);
        xu= box(:, 2);
        yl= box(:, 3);
        yu= box(:, 4);
end

% arrayfun(@(xl,xu,yl,yu) rectangle('Position',[xl yl xu-xl yu-yl],...
%     'LineWidth',0.5,'LineStyle','-',...
%     'EdgeColor',varargin{2},'FaceColor',varargin{1}),...
%     xl, xu, yl, yu) ;
arrayfun(@(xl,xu,yl,yu) patch('Faces',[1 2 3 4],...
    'Vertices',[xl yl; xu yl; xu yu; xl yu],...
    'LineWidth',0.5,'LineStyle','-',...
    'EdgeColor',varargin{2}, 'EdgeAlpha',varargin{3},...
    'FaceColor',varargin{1}, 'FaceAlpha',varargin{3}),...
    xl, xu, yl, yu);

end
