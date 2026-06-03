fluidPage(

    titlePanel("WHO Happiness Map app"),

    sidebarLayout(
        sidebarPanel(
            selectInput(
              "region",
              "Select Region: ",
              unique(combo_df$region),
              "Western Europe"
            ),
            selectInput(
              "variable",
              "Select Variable: ",
              combo_df %>% 
                dplyr::select(ladder_score:perceptions_of_corruption) %>% 
                colnames(),
              "ladder_score"
            )
        ),

        mainPanel(
            girafeOutput("map")
        )
    )
)
