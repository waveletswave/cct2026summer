function(input, output, session) {

    output$map <- renderGirafe({

      gg_map <- combo_df %>% 
        dplyr::filter(region == input$region) %>%
        dplyr::select(country_name, long, lat, group, input$variable) %>% 
        dplyr::rename(chosen_var = 5) %>% 
        ggplot(aes(x = long, 
                   y = lat,
                   group = group)) + # so that only lat and long points for each region are connected
        geom_polygon_interactive(aes(fill = chosen_var, 
                                     tooltip = sprintf("%s<br/>%s", country_name, chosen_var), 
                                     data_id = country_name), 
                                 color = "white") + # take a look at ?geom_polygon example for more info
        coord_sf() +
        theme_void() +
        labs(fill = input$variable) +
        theme(legend.position = "bottom")
      
      girafe(ggobj = gg_map)

    })

}
