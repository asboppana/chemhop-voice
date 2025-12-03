import containerQueries from "@tailwindcss/container-queries";
import tailwindcssAnimate from "tailwindcss-animate";


/** @type {import('tailwindcss').Config} */
export default {
  darkMode: ["class"],
  content: [
    "./index.html",
    './pages/**/*.{ts,tsx}',
    './components/**/*.{ts,tsx}',
    './src/**/*.{ts,tsx,js,jsx}',
    './src/main/**/*.{ts,tsx,js,jsx}',
  ],
  prefix: "",
  theme: {
  	container: {
  		center: true,
  		padding: '2rem',
  		screens: {
  			'2xl': '1400px'
  		}
  	},
  	extend: {
  		colors: {
			// Core colors
			'white': '#FFFFFF',
			'black': '#0A0A0A',
			'sedona': {
				'red': '#932610',
				'black': '#0A0A0A',
				'ai': "#214CA0",
			},
			
			// Gray scale
			'gray': {
				100: '#F9F9F9',
				150: '#F7F7F7',
				200: '#F3F3F3',
				250: '#F3F3F3',
				300: '#EDEDED',
				500: '#E8E8E6',
				700: '#C7C6C3',
				750: '#C7C6C2',
				800: '#9E9D9A',
				900: '#252525',
				1000: '#9E9D9A',
				1250: '#EDEDED',
				1500: '#2B2B2B',
				1750: '#252525',
				2000: '#1C1C1C',
			},
			'warm-gray': {
				500: '#9D9284',
				1000: "#5E5648",
			},
			'cool-gray': {
				200: '#DDDDDD',
				400: "#C1BFC9",
				600: "#78758A",
				800: "#4E4B63",
			},
			'purple': {
				200: '#F4E2F7',
				800: "#A08FC0",
				1000: "#8C6FC1",
			},
			// Biomarker colors
			'biomarker': {
				'green': '#4EAC72',
				'red': {
					'on-black': '#CB5A43',
					'on-white': '#B43F27',
				},
				'purple': {
					'on-black': '#A08FC0',
					'on-white': '#8C6FC1',
				},
				'orange': '#E39542',
			},

			'tags': {
				'exercise': '#3E8A9340',
				'nutrition': '#62344140',
				'supplementation': '#4A398D40',
				'clinical': '#26518640',
				'mindfulness': '#8F471340',
				'sleep': '#232A3640',
				'environment': '#CAE1C980',
				'social': '#ADDDCF80',
				'lifestyle': '#ADDDCF80',
			},

			'status-tags': {
				'draft-bg': '#E16BF826',
				'draft': '#433845',
				'active-bg': '#77F86B33',
				'active': '#2D5529',
				'quick-win': '#CFF4D4',
			},
			
			// HRZ (Heart Rate Zone) colors - from VO2DetailView
			'hrz': {
				1: '#DDDDDD', // lightest
				2: '#C1BFC9',
				3: '#78758A',
				4: '#4E4B63',
				5: '#0F0E15', // darkest
			},
			
			// Sleep stage colors - from EventDetailModal
			'sleep': {
				'awake': '#e2e8f0',
				'light': '#a78bfa',
				'deep': '#7c3aed',
				'rem': '#8b5cf6',
			},
			'dangerous': '#df1414',

			// Shadcn/ui colors
  			border: "hsl(var(--border))",
  			input: "hsl(var(--input))",
  			ring: "hsl(var(--ring))",
  			background: "hsl(var(--background))",
  			foreground: "hsl(var(--foreground))",
  			primary: {
  				DEFAULT: "hsl(var(--primary))",
  				foreground: "hsl(var(--primary-foreground))"
  			},
  			secondary: {
  				DEFAULT: "hsl(var(--secondary))",
  				foreground: "hsl(var(--secondary-foreground))"
  			},
  			destructive: {
  				DEFAULT: "#B43F27",
  				foreground: "#ffffff"
  			},
  			muted: {
  				DEFAULT: "hsl(var(--muted))",
  				foreground: "hsl(var(--muted-foreground))"
  			},
  			accent: {
  				DEFAULT: "hsl(var(--accent))",
  				foreground: "hsl(var(--accent-foreground))"
  			},
  			popover: {
  				DEFAULT: "hsl(var(--popover))",
  				foreground: "hsl(var(--popover-foreground))"
  			},
  			card: {
  				DEFAULT: "hsl(var(--card))",
  				foreground: "hsl(var(--card-foreground))"
  			},
  			chart: {
  				'1': 'hsl(var(--chart-1))',
  				'2': 'hsl(var(--chart-2))',
  				'3': 'hsl(var(--chart-3))',
  				'4': 'hsl(var(--chart-4))',
  				'5': 'hsl(var(--chart-5))'
  			},
  		},
  		borderRadius: {
  			lg: "var(--radius)",
  			md: "calc(var(--radius) - 2px)",
  			sm: "calc(var(--radius) - 4px)"
  		},
  		keyframes: {
  			"accordion-down": {
  				from: { height: 0 },
  				to: { height: "var(--radix-accordion-content-height)" },
  			},
  			"accordion-up": {
  				from: { height: "var(--radix-accordion-content-height)" },
  				to: { height: 0 },
  			},
  			"shimmer": {
  				"0%": { transform: "translateX(-100%) skewX(-12deg)" },
  				"100%": { transform: "translateX(300%) skewX(-12deg)" }
  			}
  		},
  		animation: {
  			"accordion-down": "accordion-down 0.2s ease-out",
  			"accordion-up": "accordion-up 0.2s ease-out",
	  				"shimmer": "shimmer 2s linear infinite"
  		},
  		fontFamily: {
			'pp-neue': ['PP Neue Montreal', 'sans-serif'],
			'pp-neue-mono': ['PP Neue Montreal Mono', 'monospace'],
			'monument': ['ABC Monument Grotesk Mono', 'monospace'],
  		},
  		fontWeight: {
      'thin': '100',
      'book': '400',
      'medium': '500',
      'bold': '700',
    },
  	}
  },
  plugins: [tailwindcssAnimate, containerQueries],
  safelist: [
    'bg-biomarker-green',
    'bg-biomarker-red-on-white',
    'bg-biomarker-purple-on-white',
	'bg-biomarker-red-on-black',
	'bg-biomarker-purple-on-black',
	// Add tags colors
	'bg-tags-exercise', 'bg-tags-nutrition', 'bg-tags-supplementation', 'bg-tags-clinical', 'bg-tags-mindfulness', 'bg-tags-sleep', 'bg-tags-environment', 'bg-tags-social', 'bg-tags-lifestyle',
	// tags with opacity variants used in UI
	'bg-tags-exercise/20', 'bg-tags-nutrition/20', 'bg-tags-supplementation/20', 'bg-tags-clinical/20', 'bg-tags-mindfulness/20', 'bg-tags-sleep/20', 'bg-tags-environment/20', 'bg-tags-social/20', 'bg-tags-lifestyle/20',
	// text colors for tags
	'text-tags-exercise', 'text-tags-nutrition', 'text-tags-supplementation', 'text-tags-clinical', 'text-tags-mindfulness', 'text-tags-sleep', 'text-tags-environment', 'text-tags-social', 'text-tags-lifestyle',
	// Add status tags colors
	'bg-status-tags-draft', 'bg-status-tags-active', 'bg-status-tags-quick-win',
	// Add warm gray colors
	'bg-warm-gray-500',
	// Add HRZ colors
	'text-hrz-1', 'text-hrz-2', 'text-hrz-3', 'text-hrz-4', 'text-hrz-5',
	'bg-hrz-1', 'bg-hrz-2', 'bg-hrz-3', 'bg-hrz-4', 'bg-hrz-5',
	// Add sleep colors
	'bg-sleep-awake', 'bg-sleep-light', 'bg-sleep-deep', 'bg-sleep-rem',
  ],
}

