# data file for the Fltk User Interface Designer (fluid)
version 1.0300 
header_name {.h} 
code_name {.cxx}
class TestUI {open
} {
  Function {TestUI()} {open
  } {
    Fl_Window {} {open
      xywh {1159 481 798 262} type Double resizable visible
    } {
      Fl_Scroll {} {
        label Turbulence open
        xywh {15 30 400 144} box DOWN_BOX labelsize 10 align 17
      } {
        Fl_Spinner {} {
          label {coef. et }
          xywh {98 50 87 25} type Float minimum 0 maximum 0.005 step 1e-005 value 0.001
        }
        Fl_Spinner {} {
          label {max Et }
          xywh {98 80 87 25} type Float minimum 0 maximum 0.2 step 0.0002 value 0.1
        }
        Fl_Spinner {} {
          label {sp_et }
          xywh {292 50 87 25} type Float minimum 0 maximum 200 step 0.2 value 100
        }
        Fl_Spinner {} {
          label {sp_et_mesh }
          xywh {292 80 87 25} type Float minimum 0 maximum 20 step 0.02 value 0.5
        }
        Fl_Spinner {} {
          label {threshold }
          xywh {98 140 87 25} minimum -2000 maximum 2000 step 10 value 700
        }
        Fl_Spinner {} {
          label {w. scale }
          xywh {98 110 87 25} type Float minimum 0 maximum 10 step 0.1
        }
        Fl_Spinner {} {
          label sp_et_cri
          xywh {292 110 87 25} type Float minimum 0 maximum 1 step 0.001 value 0.5
        }
        Fl_Button {} {
          label Apply
          xywh {292 140 87 25}
        }
      }
      Fl_Scroll {} {
        label Draw open
        xywh {425 31 260 144} box DOWN_BOX labelsize 10 align 17 resizable
      } {
        Fl_Value_Slider {} {
          label {Veloc. Scale }
          xywh {551 55 120 25} type Horizontal align 4 maximum 0.1 step 0.001
        }
        Fl_Value_Slider {} {
          label {Mesh Threshold }
          xywh {551 85 120 25} type Horizontal align 4 maximum 0.1 step 0.001
        }
        Fl_Check_Button {} {
          label {Refraction } selected
          xywh {550 115 25 25} down_box DOWN_BOX align 4
        }
      }
    }
  }
} 
