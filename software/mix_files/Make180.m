function ang = Make180(ang)

    ang = angle(exp(1i * ang*pi/180)) * 180/pi;