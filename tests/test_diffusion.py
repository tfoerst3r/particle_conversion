import pytest
from particle_conversion.models.diffusion import DiffusionCoeff
from math import trunc


class Test_DiffusionCoeff:

    @pytest.fixture
    def base_input(self) -> None:
        self.temp = 296.1 #K
        self.pressure = 1.01325e5 #Pa
        self.composition = {'CO':0,'CO2':1}
        self.diffobject = DiffusionCoeff()
    
    def test_diffusion(self,base_input):
        value = self.diffobject.diffusion_mix(comp=self.composition,temp=self.temp, pressure=self.pressure)
        value = trunc(value*10**9)
        assert value == 15129
    
    def test_binarydiffusion(self,base_input):
        gases = list(self.composition.keys())
        value = self.diffobject.binarydiffusion(temp=self.temp,pressure=self.pressure,gasA=gases[0],gasB=gases[1])
        value = trunc(value*10**9)
        assert value == 15129

