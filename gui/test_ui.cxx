// generated by Fast Light User Interface Designer (fluid) version 1.0300

#include "test_ui.h"

TestUI::TestUI() {
  Fl_Double_Window* w;
  { Fl_Double_Window* o = new Fl_Double_Window(465, 262);
    w = o;
    o->user_data((void*)(this));
    { Fl_Scroll* o = new Fl_Scroll(15, 30, 400, 144, "Turbulence");
      o->box(FL_DOWN_BOX);
      o->labelsize(10);
      o->align(Fl_Align(FL_ALIGN_TOP|FL_ALIGN_INSIDE));
      { Fl_Spinner* o = new Fl_Spinner(98, 50, 87, 25, "coef. et ");
        o->type(1);
        o->minimum(0);
        o->maximum(0.005);
        o->step(1e-005);
        o->value(0.001);
      } // Fl_Spinner* o
      { Fl_Spinner* o = new Fl_Spinner(98, 80, 87, 25, "max Et ");
        o->type(1);
        o->minimum(0);
        o->maximum(0.2);
        o->step(0.0002);
        o->value(0.1);
      } // Fl_Spinner* o
      { Fl_Spinner* o = new Fl_Spinner(292, 50, 87, 25, "sp_et ");
        o->type(1);
        o->minimum(0);
        o->maximum(200);
        o->step(0.2);
        o->value(100);
      } // Fl_Spinner* o
      { Fl_Spinner* o = new Fl_Spinner(292, 80, 87, 25, "sp_et_mesh ");
        o->type(1);
        o->minimum(0);
        o->maximum(20);
        o->step(0.02);
        o->value(0.5);
      } // Fl_Spinner* o
      { Fl_Spinner* o = new Fl_Spinner(98, 140, 87, 25, "threshold ");
        o->minimum(-2000);
        o->maximum(2000);
        o->step(10);
        o->value(700);
      } // Fl_Spinner* o
      { Fl_Spinner* o = new Fl_Spinner(98, 110, 87, 25, "w. scale ");
        o->type(1);
        o->minimum(0);
        o->maximum(10);
        o->step(0.1);
      } // Fl_Spinner* o
      { Fl_Spinner* o = new Fl_Spinner(292, 110, 87, 25, "sp_et_cri");
        o->type(1);
        o->minimum(0);
        o->maximum(1);
        o->step(0.001);
        o->value(0.5);
      } // Fl_Spinner* o
      { new Fl_Button(292, 140, 87, 25, "Apply");
      } // Fl_Button* o
      o->end();
      Fl_Group::current()->resizable(o);
    } // Fl_Scroll* o
    o->end();
  } // Fl_Double_Window* o
}
