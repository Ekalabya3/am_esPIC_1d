        figure; % Create a new figure
        scatter(xp1, vp1, 5, 'red', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5); % Smaller, semi-transparent red points
        hold on
        scatter(xp2, vp2, 5, 'blue', 'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeAlpha', 0.5); % Smaller, semi-transparent blue points
        xlabel('$x$ (Position)', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('$v$ (Velocity)', 'Interpreter', 'latex', 'FontSize', 12);
        pause(0.00001)
        hold off

        % Construct the filename
        filename = sprintf('scatter_plot.png');

        % Save the figure as a PNG file
        saveas(gcf, filename);

        % Close the figure to save memory
        close(gcf);

        figure;
        histogram(vp1, 50, 'Normalization', 'count', 'FaceColor', 'r', 'EdgeColor', 'r'); % Red color for vp1
        hold on
        histogram(vp2, 50, 'Normalization', 'count', 'FaceColor', 'b', 'EdgeColor', 'b'); % Blue color for vp2
        % Add axis labels with LaTeX-style formatting
        xlabel('$v$ (Velocity)', 'Interpreter', 'latex', 'FontSize', 12);
        ylabel('Number of particles', 'Interpreter', 'latex', 'FontSize', 12);
        % Optional: Add title
        % title('Particle Velocity Distribution', 'Interpreter', 'latex', 'FontSize', 14);
        pause(0.00001);
        hold off

        % Construct the filename
        filename = sprintf('histogram_plot.png');

        % Save the figure as a PNG file
        saveas(gcf, filename);

        % Close the figure to save memory
        close(gcf);