function [time,data] = acquire(shotnum)

% Acquire the data from MDS plus

    node_string = ['\i_p            '
        '\b_p0_0_t       ' % i = 2
        '\b_p0_45_t      '
        '\b_p0_90_t      '
        '\b_p0_135_t     '
        '\b_p0_180_t     '
        '\b_p0_225_t     '
        '\b_p0_270_t     '
        '\b_p0_315_t     '
        '\b_p15_0_t      ' % i = 10
        '\b_p15_45_t     '
        '\b_p15_90_t     '
        '\b_p15_135_t    '
        '\b_p15_180_t    '
        %     '\b_p15_225_t    '
        '\b_p15_270_t    '
        '\b_p15_315_t    '
        '\b_p30_0_t      ' % i = 17
        '\b_p30_45_t     '
        '\b_p30_90_t     '
        '\b_p30_135_t    '
        '\b_p30_180_t    '
        '\b_p30_225_t    '
        '\b_p30_270_t    '
        '\b_p30_315_t    '
        '\b_p45_0_t      ' % i = 25
        '\b_p45_45_t     '
        '\b_p45_90_t     '
        '\b_p45_135_t    '
        '\b_p45_180_t    '
        '\b_p45_225_t    '
        '\b_p45_270_t    '
        '\b_p45_315_t    '
        '\b_p5_0_t       ' % i = 33
        '\b_p5_90_t      '
        '\b_p5_180_t     '
        '\b_p5_270_t     '
        '\b_p10_0_t      ' % i = 37
        '\b_p10_90_t     '
        %         '\b_p10_180_t    '
        '\b_p10_270_t    '
        '\b_p20_0_t      ' % i = 40
        '\b_p20_90_t     '
        '\b_p20_180_t    '
        '\b_p20_270_t    '
        '\b_p25_0_t      ' % i = 44
        '\b_p25_90_t     '
        '\b_p25_180_t    '
        '\b_p25_270_t    '
        '\b_p35_0_t      ' % i = 48
        '\b_p35_90_t     '
        '\b_p35_180_t    '
        '\b_p35_270_t    '
        '\b_p40_0_t      ' % i = 52
        '\b_p40_90_t     '
        '\b_p40_180_t    '
        '\b_p40_270_t    '
        '\m_0_p0         ' % i = 56
        '\m_0_p5         '
        '\m_0_p10        '
        '\m_0_p15        '
        '\m_0_p10        '
        '\m_0_p25        '
        '\m_0_p30        '
        '\m_0_p35        '
        '\m_0_p40        '
        '\m_0_p45        '
        '\m_1_p0         ' % i = 66
        '\m_1_p5         '
        '\m_1_p10        '
        '\m_1_p15        '
        '\m_1_p10        '
        '\m_1_p25        '
        '\m_1_p30        '
        '\m_1_p35        '
        '\m_1_p40        '
        '\m_1_p45        '
        '\phi_m_1_p0     ' % i = 76
        '\phi_m_1_p5     '
        '\phi_m_1_p10    '
        '\phi_m_1_p15    '
        '\phi_m_1_p10    '
        '\phi_m_1_p25    '
        '\phi_m_1_p30    '
        '\phi_m_1_p35    '
        '\phi_m_1_p40    '
        '\phi_m_1_p45    '
        %     '\w_cam_mon      ' % i = 86
        '\y_p0           ' % i = 86
        '\y_p5           '
        '\y_p10          '
        '\y_p15          '
        '\y_p20          '
        '\y_p25          '
        '\y_p30          '
        '\y_p35          '
        '\y_p40          '
        '\y_p45          '
        '\m_1_p0_norm    ' % i = 96
        '\m_1_p5_norm    '
        '\m_1_p10_norm   '
        '\m_1_p15_norm   '
        '\m_1_p10_norm   '
        '\m_1_p25_norm   '
        '\m_1_p30_norm   '
        '\m_1_p35_norm   '
        '\m_1_p40_norm   '
        '\m_1_p45_norm   '
        '\m1:mean        ' % i = 105
        '\m1:mean:error  '
        ];
    
    tree_string = 'zaphd';
    [time,data] = mds_read_object(shotnum, node_string, tree_string);