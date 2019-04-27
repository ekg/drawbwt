#include <iostream>
#include <string>
#include <vector>
#include <tuple>
#include <utility>
#include <cairo/cairo.h>
#include <cairo/cairo-svg.h>
#include <cairo/cairo-pdf.h>
#include "args.hxx"

int main(int argc, char** argv) {
    args::ArgumentParser parser("drawbwt: BWT renderer");
    args::HelpFlag help(parser, "help", "display this help menu", {'h', "help"});
    args::ValueFlag<std::string> output(parser, "FILE", "write the output PDF vector graphic to this file", {'o', "output-pdf"});
    args::ValueFlag<std::string> text(parser, "STR", "derive the BWT for this text", {'t', "text"});
    args::ValueFlag<double> image_width(parser, "FLOAT", "output image width", {'x', "image-width"});
    args::ValueFlag<double> image_height(parser, "FLOAT", "output image height", {'y', "image-height"});
    //args::Flag keep_temp_files(parser, "", "keep intermediate files generated during graph induction", {'k', "keep-temp"});
    args::Flag debug(parser, "debug", "enable debugging", {'d', "debug"});
    try {
        parser.ParseCLI(argc, argv);
    } catch (args::Help) {
        std::cout << parser;
        return 0;
    } catch (args::ParseError e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }
    if (argc==1) {
        std::cout << parser;
        return 1;
    }

    std::vector<std::pair<std::string, uint64_t> > v;
    std::string t = args::get(text);
    for (uint64_t i = 0; i <= t.size(); ++i) {
        std::string q = t.substr(i, t.size()-i) + "$" + t.substr(0, i);
        v.push_back(make_pair(q, i));
    }
    std::sort(v.begin(), v.end());

    cairo_surface_t *surface = nullptr;
	cairo_t *cr = nullptr;

    surface = cairo_pdf_surface_create(args::get(output).c_str(), args::get(image_width), args::get(image_height));
    cr = cairo_create(surface);

    cairo_text_extents_t te;

    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 18);
    
    //cairo_text_extents(cr, "BWM", &te);
    cairo_move_to(cr, 3, 25);
    cairo_show_text(cr, "BWM");

    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    double contents_start = 45;
    double vertical_step = 16;
    double bwm_width = 0;
    double y = contents_start;
    for (auto& s : v) {
        const std::string& seq = s.first;
        cairo_text_extents(cr, seq.c_str(), &te);
        bwm_width = std::max(bwm_width, te.width);
        cairo_move_to(cr, 3, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    bwm_width += 20;

    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 18);
    
    //cairo_text_extents(cr, "SA", &te);
    cairo_move_to(cr, bwm_width, 25);
    cairo_show_text(cr, "SA");

    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);

    double sa_width = 0;
    y = contents_start;
    for (auto& s : v) {
        const std::string& seq = std::to_string(s.second);
        cairo_text_extents(cr, seq.c_str(), &te);
        sa_width = std::max(sa_width, te.width);
        cairo_move_to(cr, bwm_width, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    // now the fun part: the F vector, the L vector (the BWT) and the LF mapping!
    
    
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    
    return(0);
}
