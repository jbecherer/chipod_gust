% returns a mask (badMotion) which is 1 when the chipod is
% presumably moving into its own wake (in x-y plane)

function [badMotion] = make_flushing_mask(motion, mask_spd, vel, do_plot)


    disp(['Generating orientation mask. Make sure compass is calibrated properly!'])

    if mask_spd == 'm'
        % more sophisticated use compass
        cmp = angle(exp(1i*motion.cmp/180*pi));

        u = interp1(vel.time, vel.u, motion.time);
        v = interp1(vel.time, vel.v, motion.time);
        velx = u.*sin(cmp) + v.*cos(cmp);
        vely = u.*cos(cmp) + v.*sin(cmp);
        backangle = atan2(v, u) * 180/pi;

        % compass == 0 -> North
        % rotate so that compass == 0 -> East
        cmpew = angle(exp(1i * -(motion.cmp-90) * pi/180))*180/pi;

        relvelx = motion.a_vel_x - velx;
        relMotionAngle = Make180(motion.chiangle - backangle);
        relSensorAngle = Make180(cmpew - backangle);

        if do_plot
            CreateFigure;
            set(gcf, 'DefaultHistogramEdgeColor', 'None')
            subplot(211)
            histogram(cmpew, 'Normalization', 'pdf'); hold on;
            histogram(motion.chiangle, 'Normalization', 'pdf');
            histogram(backangle, 'Normalization', 'pdf');
            legend('rotated compass (0 = East)', ...
                   'chipod velocity angle', ...
                   'background flow angle');
            xlim([-1,1]*180)
            xlabel('degrees'); ylabel('PDF');

            subplot(212)
            histogram(relSensorAngle, 'Normalization', 'pdf'); hold on;
            histogram(relMotionAngle, 'Normalization', 'pdf');
            legend('rel sensor angle = sensor angle - flow angle', ...
                   'rel motion angle = chipod motion angle - flow angle')
            xlim([-1,1]*180)
            xlabel('degrees'); ylabel('PDF');
        end

        sensorPointsIntoFlow = abs(relSensorAngle) > 90;
        sensorMovesInFlowDir = abs(relvelx) < abs(velx);
        sensorFasterThanFlow = abs(motion.a_vel_x) > abs(velx);

        % 1. χpod pointing into background flow
        %   a. χpod moving opp to dir of flow = ALWAYS GOOD
        %   b. χpod moving in direction of flow
        %      i. slower than flow = GOOD
        %     ii. faster than flow = BAD

        % 2. χpod pointing along background flow
        %   a. χpod moving opp to dir of flow = ALWAYS BAD
        %   b. χpod moving in direction of flow
        %      i. slower than flow = BAD
        %     ii. faster than flow = GOOD

        badMotion = logical(zeros(size(motion.time)));

        % 1 b ii
        badMotion(  sensorPointsIntoFlow ...
                    & sensorMovesInFlowDir ...
                    & sensorFasterThanFlow) = 1;

        % 2 a
        badMotion(  ~sensorPointsIntoFlow ...
                    & ~sensorMovesInFlowDir) = 1;

        % 2 b i
        badMotion( ~sensorPointsIntoFlow ...
                   & sensorMovesInFlowDir ...
                   & ~sensorFasterThanFlow) = 1;

    end
