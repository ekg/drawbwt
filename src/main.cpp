#include <iostream>
#include <string>
#include <vector>
#include <map>
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

    std::vector<char> f_v, l_v;
    for (auto& s : v) {
        f_v.push_back(s.first.at(0));
        l_v.push_back(s.first.at(s.first.size()-1));
    }
    std::vector<char> enc_v;
    std::map<char, uint64_t> enc_m;
    std::vector<uint64_t> c_v;
    uint64_t j = 0;
    for (auto& c : f_v) {
        if (enc_v.empty() || enc_v.back() != c) {
            enc_m[c] = enc_v.size();
            enc_v.push_back(c);
            c_v.push_back(j);
        }
        ++j;
    }
    c_v.push_back(l_v.size());

    cairo_surface_t *surface = nullptr;
	cairo_t *cr = nullptr;

    surface = cairo_pdf_surface_create(args::get(output).c_str(), args::get(image_width), args::get(image_height));
    cr = cairo_create(surface);

    cairo_text_extents_t te;

    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 16);
    cairo_text_extents(cr, "BWM", &te);
    double bwm_width = 0;
    bwm_width = std::max(bwm_width, te.width);
    cairo_move_to(cr, 3, 25);
    cairo_show_text(cr, "BWM");

    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    double contents_start = 45;
    double vertical_step = 16;
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
    cairo_set_font_size(cr, 16);
    cairo_text_extents(cr, "SA", &te);
    double sa_width = std::max(sa_width, te.width);
    cairo_move_to(cr, bwm_width, 25);
    cairo_show_text(cr, "SA");

    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);

    y = contents_start;
    for (auto& s : v) {
        const std::string& seq = std::to_string(s.second);
        cairo_text_extents(cr, seq.c_str(), &te);
        sa_width = std::max(sa_width, te.width);
        cairo_move_to(cr, bwm_width, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    /*
    double f_width = 0;
    double f_start_x = bwm_width + sa_width + 15;
    
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 16);
    cairo_text_extents(cr, "F", &te);
    f_width = std::max(f_width, te.width);
    cairo_move_to(cr, f_start_x, 25);
    cairo_show_text(cr, "F");
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    y = contents_start;
    for (auto& c : f_v) {
        const std::string& seq = std::string(1, c);
        cairo_text_extents(cr, seq.c_str(), &te);
        f_width = std::max(f_width, te.width);
        cairo_move_to(cr, f_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    double l_width = 0;
    double l_start_x = f_start_x + 25;
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 16);
    cairo_text_extents(cr, "L", &te);
    l_width = std::max(l_width, te.width);
    cairo_move_to(cr, l_start_x, 25);
    cairo_show_text(cr, "L");
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    
    y = contents_start;
    for (auto& c : l_v) {
        const std::string& seq = std::string(1, c);
        cairo_text_extents(cr, seq.c_str(), &te);
        l_width = std::max(l_width, te.width);
        cairo_move_to(cr, l_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }
    */

    double c_width = 0;
    //double c_start_x = l_start_x + 25;
    double c_start_x = bwm_width + sa_width + 15;
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 16);
    cairo_text_extents(cr, "C", &te);
    c_width = std::max(c_width, te.width);
    cairo_move_to(cr, c_start_x, 25);
    cairo_show_text(cr, "C");
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    
    y = contents_start;

    for (auto& c : enc_v) {
        const std::string& seq = std::string(1, c) + "→" + std::to_string(c_v[enc_m[c]]);
        cairo_text_extents(cr, seq.c_str(), &te);
        c_width = std::max(c_width, te.width);
        cairo_move_to(cr, c_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    // now LF
    // LF (i) = C[BWT[i]] + rank BWT[i] (BWT, i)
    auto lf_f = [&](uint64_t i) {
        char c = l_v[i];
        uint64_t rank_c = 0;
        for (uint64_t j = 0; j <= i; ++j) {
            if (l_v[j] == c) ++rank_c;
        }
        //std::cerr << c << " " << enc_m[c] << " from " << i << " C[c]=" << c_v[enc_m[c]] << " rank=" << rank_c << " LF = " << c_v[enc_m[c]] + rank_c << std::endl;
        return c_v[enc_m[c]] + rank_c - 1;
    };

    // draw the LF-mapping
    double lf_mapping_start_x = c_start_x + c_width + 15;
    double lf_mapping_width = 0;
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 16);
    cairo_text_extents(cr, "F", &te);
    lf_mapping_width = std::max(lf_mapping_width, te.width);
    cairo_move_to(cr, lf_mapping_start_x, 25);
    cairo_show_text(cr, "F");

    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    y = contents_start;
    for (auto& c : f_v) {
        const std::string& seq = std::string(1, c);
        cairo_text_extents(cr, seq.c_str(), &te);
        lf_mapping_width = std::max(lf_mapping_width, te.width);
        cairo_move_to(cr, lf_mapping_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }
    lf_mapping_start_x += lf_mapping_width + 60;

    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 16);
    cairo_text_extents(cr, "L", &te);
    lf_mapping_width = std::max(lf_mapping_width, te.width);
    cairo_move_to(cr, lf_mapping_start_x, 25);
    cairo_show_text(cr, "L");
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    y = contents_start;
    lf_mapping_width = 0;
    j = 0;
    for (auto& c : l_v) {
        const std::string& seq = std::string(1, c);
        cairo_text_extents(cr, seq.c_str(), &te);
        lf_mapping_width = std::max(lf_mapping_width, te.width);
        cairo_move_to(cr, lf_mapping_start_x, y);
        cairo_show_text(cr, seq.c_str());
        // draw the line
        double lf_m = lf_f(j++);
        //cairo_set_source_rgb(cr, 0, 0, 0);
        cairo_set_line_width(cr, 0.5);
        cairo_move_to(cr, lf_mapping_start_x - 2, y - 3.5);
        cairo_line_to(cr, lf_mapping_start_x - 57, contents_start + (vertical_step * lf_m) - 3.5);
        cairo_stroke(cr);
        y += vertical_step;
    }

    double lf_width = 0;
    double lf_start_x = lf_mapping_start_x + lf_mapping_width + 15;
    cairo_select_font_face(cr, "Arial",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_BOLD);
    cairo_set_font_size(cr, 16);
    cairo_text_extents(cr, "LF(i)", &te);
    lf_width = std::max(lf_width, te.width);
    cairo_move_to(cr, lf_start_x, 25);
    cairo_show_text(cr, "LF(i)");
    cairo_select_font_face(cr, "Consolas",
                           CAIRO_FONT_SLANT_NORMAL, CAIRO_FONT_WEIGHT_NORMAL);
    cairo_set_font_size(cr, 12);
    
    y = contents_start;
    uint64_t k = 0;
    for (auto& c : l_v) {
        const std::string& seq = std::to_string(lf_f(k++));
        cairo_text_extents(cr, seq.c_str(), &te);
        lf_width = std::max(lf_width, te.width);
        cairo_move_to(cr, lf_start_x, y);
        cairo_show_text(cr, seq.c_str());
        y += vertical_step;
    }

    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    
    return(0);
}
