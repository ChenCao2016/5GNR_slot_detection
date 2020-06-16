
class unit:

    freq_base = 480000
    time_base = 64

    symbol_per_slot = 14

    def __init__(self, numerology, Ts):

        if Ts%unit.freq_base != 0:
            raise Exception("sampling rate is not a ratio of 480kHz")

        TsUnit = int(Ts/unit.freq_base)

        self.symbol_time = 2 ** (5-numerology) / unit.freq_base
        self.symbol_sample = int(2 ** (5-numerology) * TsUnit)
        self.extend_cp_time = 512 * 2 ** (-6-numerology)/unit.freq_base
        self.extend_cp_sample = int(512 * 2 ** (-6-numerology) * TsUnit)
        self.normal_cp_time = 144 * 2 ** (-6-numerology)/unit.freq_base 
        self.normal_cp_sample = int(144 * 2 ** (-6-numerology) * TsUnit)
        self.normal_cp0_sample = self.normal_cp_time + 0.25/unit.freq_base
        self.normal_cp0_sample = self.normal_cp_sample + (TsUnit >> 2)

        self.slot_time = 1e-3/(numerology+1)
        self.slot_per_subframe = 1 << numerology
        self.slot_per_frame = self.slot_per_subframe * 10
        self.scs = 15000 << numerology
