function maybe_symlog(ax, axis, linthresh)
% sets symlog transform on axis if 30% of data is of opposite sign to the
% remaining 70%
%
% Inputs:
%     ax : Axes handle
%     axis : 'x' or 'y'
%     linthresh : log10(threshold within which axis is linear scaled).

    children = get(ax, 'Children');

    is_labeled = 0;

    for ccc=1:length(children)
        child = children(ccc);
        try
            if axis == 'y'
                cc = child.YData;
            else
                cc = child.XData;
            end
        catch ME
            continue
        end

        if length(cc) < 3, continue; end

        sgn = sign(nanmedian(cc));

        if sgn >= 0
            fraction = nansum(cc < 0) / length(cc);
        else
            fraction = nansum(cc > 0) / length(cc);
        end

        if fraction > 0.3
            symlog(ax, axis, linthresh);
            if ~is_labeled
                htitle = get(ax, 'Title');
                htitle.String = [htitle.String ' symmetric log axis!'];
                is_labeled = 1;
            end
        end
    end
end